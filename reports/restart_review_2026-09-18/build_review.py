#!/usr/bin/env python3
"""Build retrospective figures from saved results only; never launch MCFOST."""
from __future__ import annotations

import base64
import csv
import hashlib
import html
import json
import os
import re
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(HERE / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
import numpy as np
from astropy.table import Table
from paper_render import build_paper_pdf

REF = ROOT / 'reference'
OBS = REF / 'observations'
PRED = REF / 'predictions'
BLUE, GREEN, ORANGE, GREY = '#2369a2', '#13725b', '#cf6a32', '#888888'
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
                     'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'figure.facecolor': 'white'})
SCORE = 'fixed_normalization_region_balanced_log_rms_dex'
SCORE_PREFIX = 'strict_continuum_cavity_mass_q_source128k_uniform_v2_v1'
INPUTS: set[Path] = set()


def load_json(p: Path):
    INPUTS.add(p)
    return json.loads(p.read_text())


def table(p: Path):
    INPUTS.add(p)
    return Table.read(p, format='ascii.ecsv')


def save(fig, name, pdf):
    fig.savefig(HERE / f'{name}.png', bbox_inches='tight')
    fig.savefig(HERE / f'{name}.pdf', bbox_inches='tight')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def clean_axis(ax):
    ax.grid(alpha=.18)
    ax.set_axisbelow(True)


def coverage(audit, pdf):
    stages = []
    stage_labels = {
        'inclination_nearir_test': 'Inclination test',
        'cavity_opening_nearir_test': 'Cavity-opening test',
        'envelope_mass_nearir_test': 'Envelope-mass test',
        'dust_grid_screen_summary': 'Initial dust grid',
        'dust_refinement_screen_summary': 'Dust-grid refinement',
        'feature_grid_screen_summary': 'Feature-inclusive grid',
        'literature_screen_summary': 'Geometry–heating grid',
        'literature_boundary_summary': 'Geometry–heating extension',
        'strict_continuum_screen': 'M2: initial continuum grid',
        'strict_continuum_refinement_v1': 'M2: continuum refinement',
        'strict_continuum_3d_q_v1': 'M2: mass–size–q grid',
        'strict_continuum_cavity_mass_q_v1': 'M2: mass–cavity–q grid',
    }
    for s in audit['history_superseded']:
        if s['stage'] in {'full_R100_finalists_summary', 'full_R100_aperture_finalists_summary', 'high_photon_feature_sed_summary'}:
            continue
        stages.append((stage_labels.get(s['stage'], s['stage']), s['axes'], GREY))
    for s in audit['method2_stages']:
        axes = dict(s['axes'])
        axes.update(inclination_deg=[70], target_total_lsun=[1.9])
        if 'requested_half_opening_deg' not in axes:
            axes['cavity_half_opening_deg'] = [20]
        axes.setdefault('envelope_size_exponent', [3.5])
        stages.append((stage_labels.get(s['stage'], s['stage']), axes, BLUE))
    stages.append(('Ice-grain shape and size', {'mass_factor':[2.25], 'envelope_amax_um':[.4, 1, 2],
                  'envelope_size_exponent':[2.75], 'cavity_half_opening_deg':[17.5],
                  'inclination_deg':[50,60,70], 'target_total_lsun':[1.9]}, GREEN))
    specs = [('Envelope mass factor', ['mass_factor','envelope_mass_factor'], False),
             (r'Envelope $a_{max}$ ($\mu$m)', ['amax_um','envelope_amax_um'], True),
             (r'Grain exponent $q$', ['envelope_size_exponent'], False),
             ('Cavity half-opening (deg)', ['requested_half_opening_deg','cavity_half_opening_deg'], False),
             ('Inclination (deg from pole)', ['inclination_deg'], False),
             ('Source luminosity target\n'+r'($L_\odot$)', ['target_total_lsun'], False)]
    fig, axs = plt.subplots(1, 6, figsize=(17, 8.4), sharey=True)
    for ax, (label, keys, log) in zip(axs, specs):
        for i, (_, axes, color) in enumerate(stages):
            vals = next((axes[k] for k in keys if k in axes), [])
            if vals:
                ax.plot([min(vals),max(vals)], [i,i], color=color, lw=1.1, alpha=.6)
                ax.scatter(vals, np.full(len(vals), i), c=color, s=28, zorder=3)
        if log:
            ax.set_xscale('log'); ax.set_xticks([.1,.4,1,3], ['0.1','0.4','1','3'])
        ax.set_xlabel(label); clean_axis(ax)
        ax.tick_params(axis='y', length=0)
    axs[0].set_yticks(range(len(stages)), [s[0] for s in stages])
    axs[0].invert_yaxis()
    fig.suptitle('Sampled parameter space', fontsize=14)
    fig.text(.29,.025,'Grey: earlier estimators    Blue: Method-2 (M2) continuum grids    Green: ice-grain calculations\n'
             'Symbols denote sampled values. Lines indicate the sampled range, not continuous coverage. Fixed values are included where documented.', fontsize=9)
    fig.subplots_adjust(left=.23, right=.98, top=.91, bottom=.15, wspace=.28)
    save(fig, '01_parameter_coverage', pdf)


def landscape(scores, pdf):
    p = scores[np.asarray(scores['scenario']) == 'primary']
    masses = sorted(set(p['mass_factor']))
    qs = sorted(set(p['envelope_size_exponent']))
    angles = sorted(set(p['requested_half_opening_deg']))
    vmin, vmax = float(min(p[SCORE])), float(max(p[SCORE]))
    fig, axs = plt.subplots(1,3,figsize=(13.5,5.2), sharey=True)
    for ax,mass in zip(axs,masses):
        z=np.full((len(qs),len(angles)),np.nan)
        near=np.zeros_like(z,dtype=bool)
        for r in p[p['mass_factor']==mass]:
            i=qs.index(r['envelope_size_exponent']); j=angles.index(r['requested_half_opening_deg'])
            z[i,j]=r[SCORE]; near[i,j]=r['unresolved_with_scenario_rank1']
        im=ax.imshow(z,origin='lower',aspect='auto',cmap='viridis_r',vmin=vmin,vmax=vmax)
        for i in range(len(qs)):
            for j in range(len(angles)):
                ax.text(j,i,f'{z[i,j]:.3f}',ha='center',va='center',color='black' if z[i,j] < .15 else 'white',fontsize=9)
                if near[i,j]:ax.add_patch(Rectangle((j-.46,i-.45),.92,.9,fill=False,ec='#df7027',lw=2.4))
        ax.set_title(f'Envelope mass factor {mass:g}')
        ax.set_xticks(range(len(angles)),[f'{x:g}' for x in angles]);ax.set_yticks(range(len(qs)),[f'{x:g}' for x in qs])
        ax.set_xlabel('Requested cavity half-opening (deg)')
    axs[0].set_ylabel(r'Grain-size exponent $q$')
    cb=fig.colorbar(im,cax=fig.add_axes([.9,.3,.015,.5]));cb.set_label('Region-balanced log RMS (dex)')
    fig.suptitle('Continuum score distribution: mass–cavity–grain-exponent grid',fontsize=14)
    fig.text(.08,.06,'Orange boxes: 8 models within the operational 0.01086-dex tolerance (not confidence intervals).\n'
             'Nominal best: 2.25 / 17.5° / q=2.75, score 0.10931 dex. Cavity 20° differs by only 0.00045 dex.\n'
             'Fixed: i=70°, amax=0.4 µm, source target 1.9 Lsun, original Mie grains with 5% ice mantle by volume.',fontsize=10)
    fig.subplots_adjust(left=.07,right=.86,top=.81,bottom=.29,wspace=.15)
    save(fig,'02_final_continuum_landscape',pdf)


def obs_grey(ax,native):
    for segment in dict.fromkeys(str(x) for x in native['segment']):
        t=native[np.asarray(native['segment'])==segment]
        w=np.asarray(t['wavelength_um']); f=np.asarray(t['lambda_f_lambda_w_m2'])
        good=np.isfinite(f)&(f>0)&np.isfinite(w)
        ax.scatter(w[good],f[good],color='#969696',s=4,alpha=.42,rasterized=True)


def continuum(pred,pdf):
    native=table(OBS/'continuum_sed_R100.ecsv')
    fig,axs=plt.subplots(2,1,figsize=(12.8,8),sharex=True,gridspec_kw={'height_ratios':[2.2,1]})
    obs_grey(axs[0],native)
    colors=[BLUE,ORANGE,'#9254a1']
    for mi,c,label in zip([1,2,26],colors,['Nominal: m2.25 / q2.75 / cavity17.5°','Near-tie: m2.25 / q2.75 / cavity20°','Second family: m2.50 / q3.25 / cavity17.5°']):
        t=pred[pred['model_index']==mi];t.sort('wavelength_um')
        w=np.asarray(t['wavelength_um']);f=np.asarray(t['model_aperture_flux_147pc_w_m2']);o=np.asarray(t['observed_flux_geometric_w_m2'])
        axs[0].plot(w,f,'o-',color=c,lw=1.8,ms=4,label=label)
        axs[1].plot(w,np.log10(f/o),'o-',color=c,lw=1.6,ms=4)
        if mi==1:axs[0].errorbar(w,o,yerr=t['observed_error_w_m2'],fmt='ko',ms=5,capsize=2,label='Nine scored JWST anchors')
    for ax in axs:
        ax.set_xscale('log');clean_axis(ax)
        for lo,hi in [(2.6,3.85),(4.1,5.25),(5.5,13.3),(14.5,16.5),(17,23)]:ax.axvspan(lo,hi,color='grey',alpha=.065)
    axs[0].set_yscale('log');axs[0].set_ylim(1e-16,2e-12);axs[0].set_xlim(.92,30)
    axs[0].set_ylabel(r'$\lambda F_\lambda$ (W m$^{-2}$), 1″ radius at 147 pc')
    axs[0].legend(fontsize=9,loc='upper left')
    axs[0].set_title('(a) Aperture-matched continuum comparison')
    axs[1].axhline(0,color='black',lw=.8);axs[1].axhspan(-.1,.1,color=BLUE,alpha=.06)
    axs[1].text(.02,.92,'(b)',transform=axs[1].transAxes,va='top')
    axs[1].set_ylabel('log10(model / JWST), dex');axs[1].set_xlabel('Wavelength (µm)')
    axs[1].set_xticks([1,1.5,2,3,5,10,20,30],[1,1.5,2,3,5,10,20,30])
    fig.text(.09,.025,'Grey points are observed R≈100 bins, including the excluded features. Connecting model segments are guides, not predictions inside those features.\n'
             'The ≈0.995-µm anchor is excluded. The score is a regional log-RMS screening statistic, not reduced chi-square or a posterior.',fontsize=9)
    fig.subplots_adjust(left=.10,right=.98,top=.94,bottom=.14,hspace=.07)
    save(fig,'03_final_continuum_spectrum',pdf)


def ice_figures(pdf):
    dust=load_json(PRED/'h2o_dust_axis_scores.json')['near_infrared_1arcsec']
    abundance=load_json(PRED/'h2o_abundance_curve.json')
    tau=load_json(PRED/'h2o_best_dust_tau.json')
    sil=table(PRED/'silicate_demo_anchors.ecsv')
    fig,axs=plt.subplots(2,2,figsize=(13.2,9.1))
    ax=axs[0,0];matrix=np.zeros((3,3))
    for r in dust['variants']:
        if float(r['inclination_deg'])==70:
            i=[.1,.4,.8].index(float(r['dhs_vmax']));j=[.4,1,2].index(float(r['envelope_amax_um']))
            matrix[i,j]=r['slope_rms_dex']
    im=ax.imshow(matrix,origin='lower',cmap='viridis_r',aspect='auto')
    for i in range(3):
        for j in range(3):ax.text(j,i,f'{matrix[i,j]:.3f}',ha='center',va='center',color='black' if matrix[i,j]<.7 else 'white')
    ax.set_xticks(range(3),['0.4','1.0','2.0']);ax.set_yticks(range(3),['0.1','0.4','0.8'])
    ax.set_xlabel('Envelope amax (µm)');ax.set_ylabel('DHS maximum hollow-volume fraction')
    ax.set_title('(a) Near-infrared continuum residuals, i=70°')
    fig.colorbar(im,ax=ax,shrink=.8).set_label('RMS (dex)')
    ax=axs[0,1];p=abundance['points']
    ax.plot([r['h2o_mass_fraction'] for r in p],[r['unweighted_rms_tau_residual'] for r in p],'o-',color=ORANGE)
    ax.set_xlabel('H₂O species mass fraction');ax.set_ylabel('Nine-anchor RMS in optical depth')
    ax.set_title('(b) H₂O abundance dependence, vmax=0.8')
    ax.annotate('Best sampled x=0.04\nRMS τ=0.192',xy=(.04,.1916),xytext=(.022,.3),arrowprops={'arrowstyle':'->'},fontsize=9)
    clean_axis(ax)
    ax=axs[1,0];w=np.asarray(tau['wavelength_um']);scored=np.asarray(tau['scored'])
    old=next(r for r in p if r['h2o_mass_fraction']==.04)
    dense_path=OBS/'h2o_v3/stellar_center_h2o_dense_R2400_v3.csv'
    INPUTS.add(dense_path)
    with dense_path.open() as stream:dense=list(csv.DictReader(stream))
    ax.scatter([float(r['wavelength_um']) for r in dense],
               [float(r['canonical_sparse_huber_optical_depth']) for r in dense],
               s=3,c='#999999',alpha=.5,rasterized=True,label='Full observed profile')
    ax.plot(w,tau['observed_tau'],'ko-',label='JWST 0.35″ radius',ms=4)
    ax.plot(w,old['tau_profile'],'o-',color=ORANGE,label='vmax0.8: RMS τ=0.192',ms=4)
    ax.plot(w,tau['tau_profile'],'o-',color=GREEN,label='vmax0.1: RMS τ=0.287',ms=4)
    ax.axvspan(3.08,3.34,color=ORANGE,alpha=.08)
    ax.scatter(w[~scored],np.asarray(tau['observed_tau'])[~scored],facecolors='white',edgecolors='black',zorder=5)
    ax.set_xlabel('Wavelength (µm)');ax.set_ylabel('Optical depth τ');ax.set_title('(c) H₂O-band profiles, x=0.04')
    ax.legend(fontsize=8);clean_axis(ax)
    ax=axs[1,1]
    ax.plot(sil['wavelength_um'],sil['jwst_observed_flux_w_m2'],'ko-',label='JWST 1″ radius')
    ax.plot(sil['wavelength_um'],sil['best_model_flux_147pc_w_m2'],'P-',color=GREEN,label='vmax0.1 variant')
    ax.set_yscale('log');ax.set_xlabel('Wavelength (µm)');ax.set_ylabel(r'$\lambda F_\lambda$ (W m$^{-2}$)')
    ax.set_title('(d) Silicate-band predictions (not fitted)')
    ax.text(.04,.61,'At 9.7 µm: model/JWST = 0.092\nFour image-only Method-2 calculations',transform=ax.transAxes,fontsize=9)
    ax.legend(fontsize=8);clean_axis(ax)
    fig.suptitle('Continuum and absorption-band response to dust properties',fontsize=14)
    fig.text(.08,.018,'Changed-dust results reuse the old continuum temperature: these are diagnostic screens, not a self-consistent global fit.\n'
             'Original continuum dust already contained 5% ice by mantle volume; the later H₂O fraction is a species mass fraction.',fontsize=9)
    fig.subplots_adjust(left=.08,right=.96,top=.9,bottom=.13,wspace=.3,hspace=.37)
    save(fig,'04_ice_dust_tradeoffs',pdf)

    native=table(OBS/'tmc1a_sed_unstitched.ecsv')
    pred=table(PRED/'continuum_predictions.ecsv')
    base=pred[pred['model_index']==1];base.sort('wavelength_um')
    best=next(v for v in dust['variants'] if v['variant_token']=='v02_vmax0p1_amax0p4')
    a=list(best['anchors'].values())
    bw=np.array([r['wavelength_um'] for r in a]);bf=np.array([r['model_flux_147pc_w_m2'] for r in a])
    bw=np.r_[bw,np.asarray(sil['wavelength_um'])];bf=np.r_[bf,np.asarray(sil['best_model_flux_147pc_w_m2'])]
    b=table(PRED/'h2o_band_1arcsec_best.ecsv');b.sort('wavelength_um')
    all_wave=np.r_[bw,np.asarray(b['wavelength_um'])]
    all_flux=np.r_[bf,np.asarray(b['signed_aperture_flux_147pc_w_m2'])]
    order=np.argsort(all_wave)
    fig,axs=plt.subplots(2,1,figsize=(13,7.6),gridspec_kw={'height_ratios':[2.2,1]},sharex=True)
    obs_grey(axs[0],native)
    axs[0].plot(base['wavelength_um'],base['model_aperture_flux_147pc_w_m2'],'s--',color=BLUE,ms=5,label='Continuum baseline: Mie + generic ice mantle (5% volume)')
    axs[0].plot(all_wave[order],all_flux[order],'-',color=GREEN,lw=2.1,label='Near-IR preferred H₂O variant: DHS vmax0.1, x=0.04 by mass')
    axs[0].plot(bw,bf,'o',color=GREEN,ms=5)
    axs[0].plot(base['wavelength_um'],base['observed_flux_geometric_w_m2'],'ko',ms=5,label='Original nine continuum anchors')
    axs[0].plot(sil['wavelength_um'],sil['jwst_observed_flux_w_m2'],'kP',ms=7,label='Unscored silicate checks')
    axs[0].set_yscale('log');axs[0].set_ylim(1e-16,2e-12);axs[0].set_ylabel(r'$\lambda F_\lambda$ (W m$^{-2}$)')
    axs[0].set_title('(a) Observed spectrum and conditional dust-model predictions')
    axs[0].legend(fontsize=8,loc='upper left')
    r=np.array([x['log10_model_over_observed_dex'] for x in a]);sw=np.asarray(sil['wavelength_um'])
    residual_wave=np.r_[[x['wavelength_um'] for x in a],sw]
    residual_value=np.r_[r,np.asarray(sil['log10_model_over_observed_dex'])]
    ro=np.argsort(residual_wave)
    axs[1].plot(residual_wave[ro],residual_value[ro],'-',color=GREEN)
    axs[1].plot(np.array([x['wavelength_um'] for x in a]),r,'o',color=GREEN)
    axs[1].plot(sw,sil['log10_model_over_observed_dex'],'P',color=GREEN)
    axs[1].axhline(0,color='black',lw=.8);axs[1].set_ylabel('log10(model / JWST), dex');axs[1].set_xlabel('Wavelength (µm)')
    axs[1].text(.02,.1,'(b)',transform=axs[1].transAxes)
    for ax in axs:
        ax.set_xscale('log');ax.set_xlim(.95,30);clean_axis(ax)
        ax.axvspan(8,12.5,color=GREY,alpha=.09)
    axs[1].set_xticks([1,2,3,5,10,20,30],[1,2,3,5,10,20,30])
    fig.text(.10,.02,'Aperture radius: 1″. Grey points: native JWST spectrum. Connecting lines do not predict unsampled wavelengths.\n'
             'H₂O optical-depth scores use a separate 0.35″ aperture. No updated 27.5-µm prediction exists for the DHS variant.',fontsize=9)
    fig.subplots_adjust(left=.10,right=.98,top=.94,bottom=.15,hspace=.08)
    save(fig,'05_latest_full_spectrum',pdf)


def literature_plot(pdf):
    fig,axs=plt.subplots(2,2,figsize=(12.5,8))
    ax=axs[0,0]
    values=[('Harsono et al. (2014)',55,10),('Aso et al. (2015)',65,0),('Aso et al. (2021)',53,0),('Continuum model',70,0)]
    for i,(name,v,e) in enumerate(values):
        err = np.array([[.3],[.2]]) if name == 'Aso et al. (2021)' else e
        ax.errorbar(v,i,xerr=err,fmt='o',color=GREEN if name!='Continuum model' else ORANGE,capsize=4)
    ax.set_yticks(range(4),[v[0] for v in values]);ax.invert_yaxis();ax.set_xlim(43,74)
    ax.set_xlabel('Inclination (deg; 0° = face-on)');ax.set_title('(a) Published inclination estimates')
    ax.set_ylim(3.65,-.4)
    ax.text(.02,.03,'Intervals are conditional on the respective source models.',transform=ax.transAxes,fontsize=8)
    clean_axis(ax)
    ax=axs[0,1]
    for i,(label,v) in enumerate([('Legacy comparison',147),('Recent adopted value',141.8)]):ax.plot(v,i,'o',color=ORANGE if i==0 else GREEN)
    ax.set_yticks([0,1],['Legacy comparison','Adopted in Aso grain study']);ax.set_xlim(138,151);ax.set_xlabel('Distance (pc)')
    ax.set_ylim(-.25,1.3)
    ax.set_title('(b) Adopted source distances')
    ax.text(.03,.05,'Inverse-square flux ratio: (147/141.8)² = 1.075',transform=ax.transAxes,fontsize=8);clean_axis(ax)
    ax=axs[1,0]
    ax.plot(1.4415,0,'o',color=BLUE);ax.plot(1.9,1,'o',color=ORANGE);ax.plot([2.3,2.7],[2,2],'-o',color=GREEN)
    ax.set_yticks([0,1,2],['Model photosphere','Model source incl. accretion','Published bolometric values']);ax.set_xlabel(r'Luminosity ($L_\odot$)');ax.set_xlim(1.2,3)
    ax.set_ylim(-.5,2.35)
    ax.set_title('(c) Luminosity definitions')
    ax.text(.02,.05,'Source heating and inferred bolometric luminosity are distinct.',transform=ax.transAxes,fontsize=8);clean_axis(ax)
    ax=axs[1,1];ax.axis('off')
    ax.text(0,.98,'(d) Scope and applicability of\n     published constraints',fontsize=11,va='top')
    notes=[('Envelope grain size','The ≥10-µm radiative-alignment argument is not a source-specific lower limit.'),
           ('Extinction','Jet AV≈17–20 and apparent silicate AV≈33 sample different sightlines and assumptions.'),
           ('Grain populations','Disk and scattering-envelope constraints probe distinct grain populations.'),
           ('Mass and cavity','Comparisons depend on opacity, gas/dust ratio, spatial extent and full/half-angle convention.')]
    y=.79
    for label,txt in notes:
        ax.text(0,y,label+':',fontsize=10,fontweight='bold',va='top')
        ax.text(0,y-.07,textwrap.fill(txt,48),fontsize=8.4,va='top');y-=.235
    fig.suptitle('Comparison with published source constraints',fontsize=14)
    fig.text(.07,.015,'References and measurement definitions are provided in the text. The displayed estimates do not constitute a combined posterior.',fontsize=9)
    fig.subplots_adjust(left=.22,right=.96,bottom=.17,top=.88,hspace=.49,wspace=.78)
    save(fig,'06_literature_comparison',pdf)


def inventory_plot(pdf):
    j=load_json(HERE/'workspace_inventory.json')
    fig,axs=plt.subplots(1,2,figsize=(12.5,5.4))
    keys=['preserve_run_inputs_products_and_provenance','preserve_compact_science_products','path_sensitive_legacy_source','regenerable_cache']
    labels=['Run inputs / selected products','Compact science results','Legacy code','Regenerable caches']
    sizes=[j['categories'][k]['bytes']/2**20 for k in keys]
    axs[0].barh(labels,sizes,color=[BLUE,GREEN,ORANGE,GREY]);axs[0].invert_yaxis();axs[0].set_xlabel('Logical file size (MiB)')
    for i,s in enumerate(sizes):axs[0].text(s+6,i,f'{s:.1f}',va='center')
    axs[0].set_xlim(0,max(sizes)*1.18);axs[0].set_title('(a) Archived data volume by category')
    counts=[j['top_level_file_extensions'].get(x,0) for x in ['.py','.sh','.md']]
    axs[1].bar(['Python','Shell / jobs','Markdown'],counts,color=[BLUE,ORANGE,GREY])
    axs[1].set_ylabel('Top-level files');axs[1].set_title('(b) Historical source-file inventory')
    for i,n in enumerate(counts):axs[1].text(i,n+5,str(n),ha='center')
    fig.text(.09,.035,f"Before housekeeping: {j['total_files']:,} files, {j['total_bytes']/2**30:.2f} GiB logical size. "
             f"{j['top_level_path_relative_source_count']} source files resolve resources relative to their own location.\n"
             'Historical code and scientific products are retained in a verified archive with their original relative paths.',fontsize=9)
    fig.subplots_adjust(left=.21,right=.97,bottom=.23,top=.88,wspace=.35)
    save(fig,'07_workspace_inventory',pdf)


def html_report():
    p=HERE/'REPORT.md'
    if not p.exists():return
    def inline(txt):
        txt=html.escape(txt)
        txt=re.sub(r'\*\*(.+?)\*\*',r'<strong>\1</strong>',txt)
        txt=re.sub(r'`([^`]+)`',r'<code>\1</code>',txt)
        return re.sub(r'\[([^]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',txt)
    parts=[];in_table=False;in_pre=False
    for line in p.read_text().splitlines():
        if line.startswith('```'):
            parts.append('</pre>' if in_pre else '<pre>');in_pre=not in_pre;continue
        if in_pre:parts.append(html.escape(line)+'\n');continue
        if line.startswith('|'):
            if not in_table:parts.append('<table>');in_table=True
            cells=[x.strip() for x in line.strip('|').split('|')]
            if all(set(x).issubset(set('-: ')) for x in cells):continue
            parts.append('<tr>'+''.join('<td>'+inline(x)+'</td>' for x in cells)+'</tr>');continue
        if in_table:parts.append('</table>');in_table=False
        m=re.fullmatch(r'!\[([^]]*)\]\(([^)]+)\)',line)
        if m:
            image_path=HERE/m.group(2)
            encoded=base64.b64encode(image_path.read_bytes()).decode()
            parts.append(f'<figure><img alt="{html.escape(m.group(1))}" src="data:image/png;base64,{encoded}"><figcaption>{html.escape(m.group(1))}</figcaption></figure>');continue
        if line.startswith('#'):
            n=len(line)-len(line.lstrip('#'));parts.append(f'<h{n}>'+inline(line[n:].strip())+f'</h{n}>')
        elif line.startswith('- '):parts.append('<p class="item">• '+inline(line[2:])+'</p>')
        elif line:parts.append('<p>'+inline(line)+'</p>')
    if in_table:parts.append('</table>')
    css='body{max-width:1000px;margin:45px auto;padding:0 32px;font:17px/1.6 Georgia,"Times New Roman",serif;color:#181818}h1{font-size:29px;line-height:1.25;text-align:center}h2{font-size:23px;margin-top:1.8em}h3{font-size:19px}p{text-align:justify}img{width:100%;height:auto}table{border-collapse:collapse;width:100%;font-size:14px;margin:20px 0}td{border-bottom:1px solid #ccc;padding:8px;vertical-align:top}tr:first-child{border-top:2px solid #444;border-bottom:1px solid #444;font-weight:600}figure{margin:30px 0 12px}figcaption{color:#444;font-size:14px}pre{white-space:pre-wrap;background:#f5f5f5;padding:15px}a{color:#285779}.item{margin-left:15px}@media print{body{font-size:11px}figure{break-inside:avoid}h2{break-after:avoid}}'
    (HERE/'REPORT.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Radiative-transfer modelling of TMC1A</title><style>'+css+'</style></head><body>'+''.join(parts)+'</body></html>')


def main():
    INPUTS.update([HERE/'REPORT.md', HERE/'build_review.py', HERE/'paper_render.py'])
    audit=load_json(HERE/'continuum_audit.json')
    scores=table(PRED/'continuum_scores.ecsv')
    pred=table(PRED/'continuum_predictions.ecsv')
    with PdfPages(HERE/'FIGURES.pdf') as pdf:
        coverage(audit,pdf);landscape(scores,pdf);continuum(pred,pdf)
        ice_figures(pdf);literature_plot(pdf);inventory_plot(pdf)
    html_report()
    build_paper_pdf(HERE)
    paths={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(INPUTS)}
    outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.png'))}
    (HERE/'figure_provenance.json').write_text(json.dumps({'date':'2026-09-18','read_only_science_inputs':True,
        'new_mcfost_runs':0,'input_sha256':paths,'figure_sha256':outputs},indent=2)+'\n')
    print('Built seven publication-style figures, FIGURES.pdf, REPORT.pdf and self-contained REPORT.html.')


if __name__=='__main__':main()
