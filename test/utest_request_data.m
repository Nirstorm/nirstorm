function data_fn = utest_request_data(data_bfn_parts)
local_fns = nst_request_files({['unittest' data_bfn_parts]}, 0);
data_fn = local_fns{1};
end