/// Derives the next positional prefix array and divergence array
/// (Durbin's a_{k+1}[i] and d_{k+1}[i]) from the current ones
///
/// # Arguments
///
/// * `pos_prefix_array` - An array slice holding the positional prefix array (a_k[i])
/// * `div_array` - An array slice holding the divergence array (d_k[i])
/// * `haplotypes` - An array slice holding the haplotypes at the current (k-th) site
fn next_pp_and_div_arrays(pos_prefix_array: &[usize], div_array: &[usize], haplotypes: &[u8], k: usize) -> (Vec<usize>, Vec<usize>) {
    let mut allele_0_pos: Vec<usize> = vec![];
    let mut allele_0_div: Vec<usize> = vec![];
    let mut allele_1_pos: Vec<usize> = vec![];
    let mut allele_1_div: Vec<usize> = vec![];

    let mut match_0_start: usize = k + 1; // if this stays k+1, no match
    let mut match_1_start: usize = k + 1;

    for (&haplo_index, &start_of_match) in pos_prefix_array.iter().zip(div_array.iter()) {
        if start_of_match > div_0_index {
            match_0_start = start_of_match;
        }
        if start_of_match > div_1_index {
            match_1_start = start_of_match;
        }
        if haplotypes[haplo_index] == 0 {
            allele_0_pos.push(haplo_index);
            allele_0_div.push(match_0_start);
            match_0_start = 0;
        } else {
            allele_1_pos.push(haplo_index);
            allele_1_div.push(match_1_start);
            match_1_start = 0;
        }
    }
    let mut next_pp_array: Vec<usize> = vec![];
    next_pp_array.extend(allele_0_pos);
    next_pp_array.extend(allele_1_pos);
    let mut next_div_array: Vec<usize> = vec![];
    next_div_array.extend(allele_0_div);
    next_div_array.extend(allele_1_div);
    (next_pp_array, next_div_array)
}


/// Return the positional prefix and divergence arrays (Durbin's a_k[i] and d_k[i])
///
/// # Arguments
///
/// * `haplotypes` - An array slice in which each row (of length n_haplo)
///                  corresponds to a single site, and contains 0s and 1s
///                  indicating allele
/// * `n_haplo` - The number of haplotypes
/// * `k_site` - The numbered site to build the array for
fn build_pp_and_div_arrays(haplotypes: &[u8], n_haplo: usize, k_site: usize) -> (Vec<usize>, Vec<usize>) {
    let mut pos_prefix_array: Vec<usize> = vec![0; n_haplo];
    for i in 0..n_haplo {
        pos_prefix_array[i] = i;
    }
    let mut div_array: Vec<usize> = vec![0; n_haplo];

    for k in 0..(k_site-1) {
        let haplotypes_at_k = &haplotypes[k*n_haplo..(k+1)*n_haplo];
        let (pos_prefix_array_tmp, div_array_tmp) = next_pp_and_div_arrays(&pos_prefix_array, &div_array, &haplotypes_at_k, k);
        pos_prefix_array = pos_prefix_array_tmp;
        div_array = div_array_tmp;
    }
    (pos_prefix_array, div_array)
}

fn main() {
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pp_and_div_arrays() {
        let n_haplo: usize = 8;
        let n_sites: usize = 6;
        let haplotypes: Vec<u8> = vec![0, 1, 1, 0, 0, 1, 1, 0,
                                       1, 1, 1, 1, 0, 0, 1, 1,
                                       0, 0, 1, 1, 0, 0, 0, 0,
                                       1, 0, 1, 1, 0, 0, 0, 1,
                                       0, 0, 1, 1, 0, 1, 0, 1,
                                       1, 1, 1, 0, 0, 0, 1, 0];
        let expected_pp_1: Vec<usize> = (0..8).collect();
        let expected_div_1: Vec<usize> = vec![0; 8];
        let (pp_1, div_1) = build_pp_and_div_arrays(&haplotypes, n_haplo, 1);
        assert_eq!(&expected_pp_1, &pp_1);
        assert_eq!(&expected_div_1, &div_1);
        let expected_pp_6: Vec<usize> = vec![4, 1, 6, 0, 5, 7, 3, 2];
        let expected_div_6: Vec<usize> = vec![5, 2, 0, 4, 5, 4, 3, 1];
        let (pp_6, div_6) = build_pp_and_div_arrays(&haplotypes, n_haplo, n_sites);
        assert_eq!(&expected_pp_6, &pp_6);
        assert_eq!(&expected_div_6, &div_6);
    }
}
