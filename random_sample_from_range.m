function val_out = random_sample_from_range(min_a, max_b, N)
    rng shuffle;    
    val_out = min_a + (max_b-min_a) .* rand(N,1);
end