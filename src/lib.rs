mod gamma_fn;
mod beta_fn;

pub use crate::gamma_fn::{log_gamma, gamma};
pub use crate::beta_fn::beta;

#[cfg(test)]
mod test_util;