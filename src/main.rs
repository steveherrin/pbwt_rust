/// Returns the next positional prefix array (Durbin's a_{k+1}) from the current one
///
/// # Arguments
///
/// * `pos_prefix_array` - An array slice holding the current positional prefix array (a_k)
/// * `haplotypes` - An array slice holding the haplotypes at the current (k-th) site
fn next_positional_prefix_array(pos_prefix_array: &[usize], haplotypes: &[u8]) -> Vec<usize> {
    let mut allele_0s: Vec<usize> = vec![];
    let mut allele_1s: Vec<usize> = vec![];

    for &i in pos_prefix_array {
        if haplotypes[i] == 0 {
            allele_0s.push(i);
        } else {
            allele_1s.push(i);
        }
    }
    let mut next_pos_prefix_array: Vec<usize> = vec![];
    next_pos_prefix_array.extend(&allele_0s);
    next_pos_prefix_array.extend(&allele_1s);
    next_pos_prefix_array
}


/// Returns the positional prefix array (Durbin's a_k)
///
/// # Arguments
///
/// * `haplotypes` - An array slice in which each row (of length n_haplo)
///                  corresponds to a single site, and contains 0s and 1s
///                  indicating allele
/// * `n_haplo` - The number of haplotypes
/// * `k_site` - The numbered site to build the array for
fn build_positional_prefix_array(haplotypes: &[u8], n_haplo: usize, k_site: usize) -> Vec<usize> {
    let mut pos_prefix_array: Vec<usize> = vec![0; n_haplo];
    for i in 0..n_haplo {
        pos_prefix_array[i] = i;
    }

    for k in 0..(k_site-1) {
        let haplotypes_at_k = &haplotypes[k*n_haplo..(k+1)*n_haplo];
        pos_prefix_array = next_positional_prefix_array(&pos_prefix_array, &haplotypes_at_k);
    }
    pos_prefix_array
}

fn main() {
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_positional_prefix_array() {
        let n_haplo: usize = 8;
        let n_sites: usize = 6;
        let haplotypes: Vec<u8> = vec![0, 1, 1, 0, 0, 1, 1, 0,
                                       1, 1, 1, 1, 0, 0, 1, 1,
                                       0, 0, 1, 1, 0, 0, 0, 0,
                                       1, 0, 1, 1, 0, 0, 0, 1,
                                       0, 0, 1, 1, 0, 1, 0, 1,
                                       1, 1, 1, 0, 0, 0, 1, 0];
        let expected_1: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6, 7];
        assert_eq!(&expected_1,
                   &build_positional_prefix_array(&haplotypes, n_haplo, 1));
        let expected_6: Vec<usize> = vec![4, 1, 6, 0, 5, 7, 3, 2];
        assert_eq!(&expected_6,
                   &build_positional_prefix_array(&haplotypes, n_haplo, n_sites));
    }
}
