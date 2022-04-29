mod gamma_fn;
mod beta_fn;
mod igamma_fn;

pub use crate::gamma_fn::{log_gamma, gamma};
pub use crate::beta_fn::beta;
pub use crate::igamma_fn::*;

#[cfg(test)]
mod test_util;