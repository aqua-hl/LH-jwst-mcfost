#!/usr/bin/env python3
"""Deterministic tests for the scale-invariant H2O continuum kernel v2."""

from __future__ import annotations

import unittest
from unittest import mock

import numpy as np

import h2o_sparse_continuum_kernel_v2 as kernel


class H2OSparseContinuumKernelV2Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.support_wave = np.asarray([2.575, 2.612, 2.650, 3.550, 3.574, 3.598])
        self.support_values = np.asarray([0.0948, 0.1043, 0.1080, 0.2578, 0.2707, 0.2725])
        self.sides = np.asarray(["blue", "blue", "blue", "red", "red", "red"])
        self.evaluation = np.linspace(2.55, 3.65, 2187)

    def fit(self, factor: float = 1.0) -> tuple[np.ndarray, dict]:
        return kernel.robust_scale_invariant_linear_continuum(
            self.evaluation,
            self.support_wave,
            factor * self.support_values,
            self.sides,
        )

    def test_scale_invariant_over_forty_orders_of_magnitude(self) -> None:
        reference, reference_record = self.fit()
        for exponent in range(-20, 21):
            factor = 10.0**exponent
            continuum, record = self.fit(factor)
            np.testing.assert_allclose(
                continuum / factor, reference, rtol=3.0e-14, atol=0.0
            )
            np.testing.assert_allclose(
                record["support_weights"],
                reference_record["support_weights"],
                rtol=3.0e-13,
                atol=3.0e-13,
            )
            self.assertAlmostEqual(
                record["dimensionless_intercept"],
                reference_record["dimensionless_intercept"],
                delta=3.0e-14,
            )
            self.assertAlmostEqual(
                record["dimensionless_slope"],
                reference_record["dimensionless_slope"],
                delta=3.0e-14,
            )

    def test_relative_convergence_and_positive_output_are_recorded(self) -> None:
        continuum, record = self.fit()
        self.assertTrue(record["converged"])
        self.assertLessEqual(
            record["final_relative_coefficient_change"],
            record["relative_coefficient_tolerance"],
        )
        self.assertGreaterEqual(record["iterations"], 1)
        self.assertLessEqual(record["iterations"], record["maximum_iterations"])
        self.assertTrue(np.all(np.isfinite(continuum) & (continuum > 0.0)))

    def test_known_slow_smooth_positive_vector_converges(self) -> None:
        # Regression supplied by the independent v3 operator audit.  The old
        # 50-iteration cap rejected this valid continuum; the unchanged
        # 1e-13 relative tolerance converges after roughly 248 iterations.
        values = np.asarray(
            [
                0.10018710,
                0.10575335,
                0.11191098,
                0.25273478,
                0.25622887,
                0.25998011,
            ]
        )
        continuum, record = kernel.robust_scale_invariant_linear_continuum(
            self.evaluation, self.support_wave, values, self.sides
        )
        self.assertTrue(record["converged"])
        self.assertGreater(record["iterations"], 50)
        self.assertLess(record["iterations"], kernel.MAX_ITERATIONS)
        self.assertTrue(np.all(np.isfinite(continuum) & (continuum > 0.0)))

    def test_second_known_slow_vector_exceeds_one_thousand_and_converges(self) -> None:
        values = np.asarray(
            [
                0.10035452,
                0.10609603,
                0.11376011,
                0.25148530,
                0.25398503,
                0.25894472,
            ]
        )
        continuum, record = kernel.robust_scale_invariant_linear_continuum(
            self.evaluation, self.support_wave, values, self.sides
        )
        self.assertTrue(record["converged"])
        self.assertGreater(record["iterations"], 1000)
        self.assertLess(record["iterations"], kernel.MAX_ITERATIONS)
        self.assertTrue(np.all(np.isfinite(continuum) & (continuum > 0.0)))

    def test_deterministic_random_smooth_positive_vectors_converge(self) -> None:
        generator = np.random.default_rng(20260726)
        normalized_x = (self.support_wave - np.mean(self.support_wave)) / np.ptp(
            self.support_wave
        )
        for _ in range(5000):
            intercept = generator.uniform(0.05, 2.0)
            slope = generator.uniform(-0.45, 0.45) * intercept
            curvature = generator.uniform(-0.08, 0.08) * intercept
            values = (
                intercept
                + slope * normalized_x
                + curvature * normalized_x**2
                + generator.normal(0.0, 0.006 * intercept, size=6)
            )
            self.assertTrue(np.all(values > 0.0))
            continuum, record = kernel.robust_scale_invariant_linear_continuum(
                self.evaluation, self.support_wave, values, self.sides
            )
            self.assertTrue(record["converged"])
            self.assertLessEqual(record["iterations"], kernel.MAX_ITERATIONS)
            self.assertTrue(np.all(np.isfinite(continuum) & (continuum > 0.0)))

    def test_iteration_cap_is_a_hard_failure(self) -> None:
        values = np.asarray(
            [
                0.10035452,
                0.10609603,
                0.11376011,
                0.25148530,
                0.25398503,
                0.25894472,
            ]
        )
        with mock.patch.object(kernel, "MAX_ITERATIONS", 1):
            with self.assertRaisesRegex(kernel.ContinuumError, "did not converge"):
                kernel.robust_scale_invariant_linear_continuum(
                    self.evaluation, self.support_wave, values, self.sides
                )

    def test_outlier_is_downweighted_without_changing_wing_roster(self) -> None:
        values = self.support_values.copy()
        values[3] *= 1.25
        _, record = kernel.robust_scale_invariant_linear_continuum(
            self.evaluation, self.support_wave, values, self.sides
        )
        weights = np.asarray(record["support_weights"])
        self.assertLess(np.min(weights), 1.0)
        self.assertEqual(len(weights), 6)

    def test_lambda_f_lambda_to_fnu_proportional(self) -> None:
        wave = np.asarray([2.0, 3.0, 4.0])
        lambda_f_lambda = np.asarray([10.0, 10.0, -2.0])
        converted = kernel.lambda_f_lambda_to_fnu_proportional(
            wave, lambda_f_lambda
        )
        np.testing.assert_array_equal(converted, [20.0, 30.0, -8.0])

    def test_invalid_support_roster_fails(self) -> None:
        with self.assertRaises(kernel.ContinuumError):
            kernel.robust_scale_invariant_linear_continuum(
                self.evaluation,
                self.support_wave,
                self.support_values,
                ["blue"] * 6,
            )
        with self.assertRaises(kernel.ContinuumError):
            kernel.robust_scale_invariant_linear_continuum(
                self.evaluation,
                self.support_wave,
                np.asarray([0.0, *self.support_values[1:]]),
                self.sides,
            )
        with self.assertRaises(kernel.ContinuumError):
            kernel.robust_scale_invariant_linear_continuum(
                [], self.support_wave, self.support_values, self.sides
            )
        invalid_wavelength = self.support_wave.copy()
        invalid_wavelength[0] = 0.0
        with self.assertRaises(kernel.ContinuumError):
            kernel.robust_scale_invariant_linear_continuum(
                self.evaluation,
                invalid_wavelength,
                self.support_values,
                self.sides,
            )

    def test_invalid_flux_conversion_fails(self) -> None:
        with self.assertRaises(kernel.ContinuumError):
            kernel.lambda_f_lambda_to_fnu_proportional([], [])
        with self.assertRaises(kernel.ContinuumError):
            kernel.lambda_f_lambda_to_fnu_proportional([2.0, 3.0], [1.0])
        with self.assertRaises(kernel.ContinuumError):
            kernel.lambda_f_lambda_to_fnu_proportional([0.0], [1.0])


if __name__ == "__main__":
    unittest.main()
