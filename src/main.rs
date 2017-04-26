fn build_prefix_array(haplotypes: &[u8], n_haplo: usize, n_sites: usize) -> Vec<usize> {
    let mut prefix_order: Vec<usize> = vec![0; n_haplo];
    for i in 0..n_haplo {
        prefix_order[i] = i;
    }

    for k in 0..(n_sites-1) {
        let mut allele_0s: Vec<usize> = vec![];
        let mut allele_1s: Vec<usize> = vec![];
        for i in &prefix_order {
            if haplotypes[i*n_sites + k] == 0 {
                allele_0s.push(*i);
            } else {
                allele_1s.push(*i);
            }
        }
        prefix_order[0..allele_0s.len()].clone_from_slice(&allele_0s);
        prefix_order[allele_0s.len()..n_haplo].clone_from_slice(&allele_1s);
    }
    prefix_order
}

fn main() {
    let n_haplo: usize = 8; // M
    let n_sites: usize = 6; // N
    let haplotypes: Vec<u8> = vec![0, 1, 0, 1, 0, 1,
                                   1, 1, 0, 0, 0, 1,
                                   1, 1, 1, 1, 1, 1,
                                   0, 1, 1, 1, 1, 0,
                                   0, 0, 0, 0, 0, 0,
                                   1, 0, 0, 0, 1, 0,
                                   1, 1, 0, 0, 0, 1,
                                   0, 1, 0, 1, 1, 0];
    assert_eq!(n_haplo * n_sites, haplotypes.len());

    let prefix_order = build_prefix_array(&haplotypes, n_haplo, n_sites);
    println!("{:?}", prefix_order);
}


