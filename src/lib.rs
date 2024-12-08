#![doc = include_str!("../README.md")]

// The following source is derived from Ryan Juckett's Damped Springs code:

/******************************************************************************
  Copyright (c) 2008-2012 Ryan Juckett
  http://www.ryanjuckett.com/

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source
     distribution.
******************************************************************************/

use num_traits::Float;

pub mod prelude {
    pub use crate::{Spring, SpringCollection, SpringConfig, SpringParams, SpringTimeStep};
}

/// Configuration options for a spring. Composed of its `angular_freq` and `damping_ratio`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpringConfig<F> {
    angular_freq: F,
    damping_ratio: F,
}

impl<F: Float> SpringConfig<F> {
    /// Construct a new spring configuration. Constrains spring parameters to non-negative values.
    pub fn new(angular_freq: F, damping_ratio: F) -> Self {
        Self {
            angular_freq: angular_freq.max(F::zero()),
            damping_ratio: damping_ratio.max(F::zero()),
        }
    }

    /// Returns the angular frequency of this spring config. Guaranteed to be at least zero.
    #[inline]
    pub fn angular_freq(&self) -> F {
        self.angular_freq
    }

    /// Returns the damping ratio of this spring config. Guaranteed to be at least zero.
    #[inline]
    pub fn damping_ratio(&self) -> F {
        self.damping_ratio
    }
}

/// Cached coefficients for a spring, based on its angular frequency and damping ratio.
///
/// Do not construct directly; instead, use [`SpringParams::from`] with a [`SpringConfig`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpringParams<F> {
    /// The spring has no angular frequency and will not move.
    Static,
    /// The spring is over-damped (damping > 1).
    OverDamped { zb: F, z1: F, z2: F },
    /// The spring is critically damped (damping = 1).
    CriticallyDamped { angular_freq: F },
    /// The spring is under-damped (damping < 1).
    UnderDamped { oz: F, a: F },
}

impl<F: Float> From<SpringConfig<F>> for SpringParams<F> {
    fn from(
        SpringConfig {
            angular_freq,
            damping_ratio,
        }: SpringConfig<F>,
    ) -> Self {
        if angular_freq < F::epsilon() {
            return Self::Static;
        }

        if damping_ratio > F::one() + F::epsilon() {
            // overdamped
            let za = -angular_freq * damping_ratio;
            let zb = angular_freq * (damping_ratio * damping_ratio - F::one()).sqrt();
            let z1 = za - zb;
            let z2 = za + zb;

            Self::OverDamped { zb, z1, z2 }
        } else if damping_ratio < F::one() - F::epsilon() {
            // under-damped
            let oz = angular_freq * damping_ratio;
            let a = angular_freq * (F::one() - damping_ratio * damping_ratio).sqrt();

            Self::UnderDamped { oz, a }
        } else {
            // critically damped
            Self::CriticallyDamped { angular_freq }
        }
    }
}

/// Cached coefficients for a spring, based on its configuration and a particular time step.
///
/// Used to efficiently update one or more springs that share the same configuration.
///
/// ## Usage
///
/// To derive a `SpringTimeStep`, you must first get a [`SpringParams`] from a [`SpringConfig`]:
///
/// ```
/// // configure your spring...
/// let config = SpringConfig {
///     angular_freq: 5,
///     damping_ratio: 0.75,
/// };
///
/// // derive its state
/// let state = SpringParams::from(config);
///
/// // then, either create your own `SpringTimeStep`...
/// let time_step = SpringTimeStep::new(state, delta_time);
/// spring.update(time_step);
///
/// // or let `Spring::update_single` do it for you
/// spring.update_single(state, delta_time);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpringTimeStep<F> {
    pp: F,
    pv: F,
    vp: F,
    vv: F,
}

impl<F: Float> Default for SpringTimeStep<F> {
    fn default() -> Self {
        Self {
            pp: F::one(),
            pv: F::zero(),
            vp: F::zero(),
            vv: F::one(),
        }
    }
}

impl<F: Float> SpringTimeStep<F> {
    /// Derive a spring time step from a particular delta time.
    pub fn new(state: impl Into<SpringParams<F>>, delta: F) -> Self {
        use SpringParams::*;
        match state.into() {
            Static => Self::default(),

            OverDamped { zb, z1, z2 } => {
                let e1 = (z1 * delta).exp();
                let e2 = (z2 * delta).exp();
                let inv_2zb = F::one() / ((F::one() + F::one()) * zb);

                let e1_2zb = e1 * inv_2zb;
                let e2_2zb = e2 * inv_2zb;

                let z1e1_2zb = z1 * e1_2zb;
                let z2e2_2zb = z2 * e2_2zb;

                Self {
                    pp: e1_2zb * z2 - z2e2_2zb + e2,
                    pv: -e1_2zb + e2_2zb,
                    vp: (z1e1_2zb - z2e2_2zb + e2) * z2,
                    vv: -z1e1_2zb + z2e2_2zb,
                }
            }

            UnderDamped { oz, a } => {
                let exp = (-oz * delta).exp();
                let cos = (a * delta).cos();
                let sin = (a * delta).sin();

                let inv_alpha = F::one() / a;

                let exp_sin = exp * sin;
                let exp_cos = exp * cos;
                let exp_ozs_alpha = exp * oz * sin * inv_alpha;

                Self {
                    pp: exp_cos + exp_ozs_alpha,
                    pv: exp_sin * inv_alpha,
                    vp: -exp_sin * a - oz * exp_ozs_alpha,
                    vv: exp_cos - exp_ozs_alpha,
                }
            }

            CriticallyDamped { angular_freq } => {
                let exp = (-angular_freq * delta).exp();
                let time_exp = delta * exp;
                let time_exp_freq = time_exp * angular_freq;

                Self {
                    pp: time_exp_freq * exp,
                    pv: time_exp,
                    vp: -angular_freq * time_exp_freq,
                    vv: -time_exp_freq * exp,
                }
            }
        }
    }

    /// Update multiple springs using this time step. Simply calls [`Spring::update`]
    /// on each passed spring.
    #[inline]
    pub fn update_many(self, springs: &mut [&mut Spring<F>]) {
        for spring in springs {
            spring.update(self);
        }
    }
}

/// An instance of a spring and its current physical properties, like its position, velocity, and target equilibrium.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spring<F> {
    pub position: F,
    pub velocity: F,
    pub equilibrium: F,
}

impl<F: Float> Default for Spring<F> {
    fn default() -> Self {
        Self {
            position: F::zero(),
            velocity: F::zero(),
            equilibrium: F::zero(),
        }
    }
}

impl<F: Float> Spring<F> {
    /// Create a new spring from a start equilibrium.
    pub fn from_equilibrium(equilibrium: F) -> Self {
        Self {
            position: F::zero(),
            velocity: F::zero(),
            equilibrium,
        }
    }

    /// Update this spring using a pre-computed [`SpringTimeStep`].
    pub fn update(&mut self, time_step: SpringTimeStep<F>) {
        let op = self.position - self.equilibrium;
        let ov = self.velocity;

        self.position = op * time_step.pp + ov * time_step.pv + self.equilibrium;
        self.velocity = op * time_step.vp + ov * time_step.vv;
    }

    /// Update this spring (and this spring only) using a [`SpringParams`] and a delta time.
    /// Will internally create a [`SpringTimeStep`] for this call.
    ///
    /// If you are updating multiple springs with the same properties,
    /// consider using [`SpringTimeStep::new`] and [`SpringTimeStep::update_many`]
    /// or create a [`SpringCollection`].
    #[inline]
    pub fn update_single(&mut self, state: SpringParams<F>, delta: F) {
        self.update(SpringTimeStep::new(state, delta));
    }
}

/// A fixed-size collection of springs that all share the same spring parameters.
/// Useful for creating springs over multiple dimensions (i.e. 2D or 3D springs).
#[derive(Debug, Clone, PartialEq)]
pub struct SpringCollection<F, const N: usize> {
    params: SpringParams<F>,
    pub springs: [Spring<F>; N],
}

impl<F: Float, const N: usize, T> From<T> for SpringCollection<F, N>
where
    T: Into<SpringParams<F>>,
{
    fn from(params: T) -> Self {
        Self {
            params: params.into(),
            springs: [Spring::default(); N],
        }
    }
}

impl<F: Float, const N: usize> SpringCollection<F, N> {
    /// Construct `N` springs, all starting at a specified `equilibrium`.
    pub fn from_equilibrium(params: impl Into<SpringParams<F>>, equilibrium: F) -> Self {
        Self {
            params: params.into(),
            springs: [Spring::from_equilibrium(equilibrium); N],
        }
    }

    /// Construct `N` springs, each with a particular equilibrium.
    pub fn from_equilibriums(params: impl Into<SpringParams<F>>, equilibriums: [F; N]) -> Self {
        Self {
            params: params.into(),
            springs: equilibriums.map(Spring::from_equilibrium),
        }
    }

    /// Update all springs over the specified delta. Constructs a new [`SpringTimeStep`]
    /// for this usage.
    #[inline]
    pub fn update(&mut self, delta: F) {
        self.update_with(SpringTimeStep::new(self.params, delta));
    }

    /// Update all springs using the specified `time_step`.
    ///
    /// **Note:** this time step need not be derived from [`SpringCollection::params`].
    /// The implementation of this method uses all constants from `time_step`.
    #[inline]
    pub fn update_with(&mut self, time_step: SpringTimeStep<F>) {
        for spring in self.springs.iter_mut() {
            spring.update(time_step);
        }
    }
}

macro_rules! impl_collection_props {
    ( $prop:ident ( $plural:ident ) => $set_prop:ident ( $set_plural:ident ) ) => {
        impl<F: Float, const N: usize> SpringCollection<F, N> {
            #[doc = concat!("The array of current spring ", stringify!($plural), ".")]
            #[inline]
            pub fn $plural(&self) -> [F; N] {
                self.springs.map(|spring| spring.$prop)
            }

            #[doc = concat!("Set a particular spring's ", stringify!($prop), ".")]
            #[inline]
            pub fn $set_prop(&mut self, index: usize, $prop: F) {
                self.springs[index].$prop = $prop;
            }

            #[doc = concat!("Set all spring ", stringify!($plural), ".")]
            #[inline]
            pub fn $set_plural(&mut self, $plural: [F; N]) {
                for (spring, $prop) in self.springs.iter_mut().zip($plural.into_iter()) {
                    spring.$prop = $prop;
                }
            }
        }
    };
}

impl_collection_props!(position (positions) => set_position (set_positions));
impl_collection_props!(velocity (velocities) => set_velocity (set_velocities));
impl_collection_props!(equilibrium (equilibriums) => set_equilibrium (set_equilibriums));
