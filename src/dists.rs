use rand::prelude::*;
//use rand::distributions::{Distribution, WeightedIndex};

/// Compute a cumulative gauge to serve as a PRNG for a truncated powerlaw $P(k)=k^{-\gamma}$.
pub fn powerlaw_cumulative_gauge(n: usize, k_min: usize, k_max: usize, gamma: f64) -> Vec<f64> {
    let dim = k_max + 1 - k_min;
    let mut k_array = vec![0; dim];
    let mut pdf = vec![0.0; dim];
    let mut cdf = vec![0.0; dim];
    k_array[0] = k_min;
    pdf[0] = n as f64 * (k_array[0] as f64).powf(-gamma);
    cdf[0] = pdf[0];
    for i in 1..dim {
        k_array[i] = k_min + i;
        pdf[i] = n as f64 * (k_array[i] as f64).powf(-gamma);
        cdf[i] = cdf[i - 1] + pdf[i];
    }
    for i in 0..dim {
        cdf[i] /= cdf[dim - 1];
    }
    cdf
}

/// Generate integer random variates using a cumulative gauge following a particular mathematical law.
pub fn random_powerlaw_variate(cdf: &[f64], k_min: usize) -> usize {
    let mut rng = rand::thread_rng();
    let trial: f64 = rng.gen();

    let mut k = 0;
    while k < cdf.len() - 1 && trial > cdf[k + 1] {
        k += 1;
    }
    k_min + k
}

/// Build a degree sequence, an array where each element is the connectivity degree of the node. This routine is
/// conceived now just for a truncated power-law function.
pub fn build_powerlaw_degree_sequence(
    n: usize,
    gamma: f64,
    k_min: usize,
    k_max: usize,
) -> Vec<usize> {
    let mut seq = vec![0; n];
    let cdf = powerlaw_cumulative_gauge(n, k_min, k_max, gamma);
    for degree in seq.iter_mut().take(n) {
        *degree = random_powerlaw_variate(&cdf, k_min);
    }
    seq
}
