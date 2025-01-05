@testset "Multi-Index Set Tests" begin
    # Test create_smolyak_miset
    n, k = 4, 3
    mi_set = create_smolyak_miset(n, k)
    @test get_n_mi(mi_set) == 35
    @test get_mi(mi_set) == [   [1, 1, 1, 1],
                                [1, 1, 1, 2],
                                [1, 1, 1, 3],
                                [1, 1, 1, 4],
                                [1, 1, 2, 1],
                                [1, 1, 2, 2],
                                [1, 1, 2, 3],
                                [1, 1, 3, 1],
                                [1, 1, 3, 2],
                                [1, 1, 4, 1],
                                [1, 2, 1, 1],
                                [1, 2, 1, 2],
                                [1, 2, 1, 3],
                                [1, 2, 2, 1],
                                [1, 2, 2, 2],
                                [1, 2, 3, 1],
                                [1, 3, 1, 1],
                                [1, 3, 1, 2],
                                [1, 3, 2, 1],
                                [1, 4, 1, 1],
                                [2, 1, 1, 1],
                                [2, 1, 1, 2],
                                [2, 1, 1, 3],
                                [2, 1, 2, 1],
                                [2, 1, 2, 2],
                                [2, 1, 3, 1],
                                [2, 2, 1, 1],
                                [2, 2, 1, 2],
                                [2, 2, 2, 1],
                                [2, 3, 1, 1],
                                [3, 1, 1, 1],
                                [3, 1, 1, 2],
                                [3, 1, 2, 1],
                                [3, 2, 1, 1],
                                [4, 1, 1, 1]];

    # Test add_mi
    mi_set_new = MISet([[1,1,1,1]]) 
    combined_mi_set = add_mi(mi_set, mi_set_new)
    @test get_n_mi(combined_mi_set) == 35
    mi_set_new = MISet([[5,1,1,1]]) 
    combined_mi_set = add_mi(mi_set, mi_set_new)
    @test get_n_mi(combined_mi_set) == 36

    # Test get_margin and get_reduced_margin
    mi_set = MISet([[1,1,1,1]]) 
    margin_set = get_margin(mi_set)
    @test get_mi(margin_set) ==  [  [1, 1, 1, 2],
                                    [1, 1, 2, 1],
                                    [1, 2, 1, 1],
                                    [2, 1, 1, 1]];
    
    reduced_margin_set = get_reduced_margin(mi_set)
    @test get_mi(reduced_margin_set) ==  [  [1, 1, 1, 2],
                                            [1, 1, 2, 1],
                                            [1, 2, 1, 1],
                                            [2, 1, 1, 1]];

    # More complicated example
    mi_set = MISet([[1, 1, 1, 1],
                    [1, 1, 1, 2],
                    [1, 1, 1, 3],
                    [1, 1, 2, 1],
                    [1, 1, 2, 2],
                    [1, 1, 2, 3],
                    [1, 1, 3, 1],
                    [1, 2, 1, 1],
                    [1, 2, 1, 2],
                    [1, 2, 2, 1],
                    [1, 2, 2, 2],
                    [1, 2, 3, 1],
                    [2, 1, 1, 1],
                    [2, 1, 1, 2],
                    [2, 1, 1, 3],
                    [2, 1, 2, 1],
                    [2, 2, 1, 1],
                    [2, 2, 1, 2],
                    [2, 2, 2, 1]]);

    rm = get_reduced_margin(mi_set);

    @test get_mi(rm) == [[1, 1, 1, 4],
                        [1, 1, 3, 2],
                        [1, 1, 4, 1],
                        [1, 2, 1, 3],
                        [1, 3, 1, 1],
                        [2, 1, 2, 2],
                        [2, 1, 3, 1],
                        [3, 1, 1, 1]];
end