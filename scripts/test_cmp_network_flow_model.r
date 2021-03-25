# test network flow model construction
source("scripts/generate_mc_mgraph_network/gen_network.r")

test_network_flow_model = function() {

    build_sing_emb_wt10_network(net_id = "test")
    message("computed network flow model")
    mct1 = scdb_mctnetwork("sing_emb_wt10")
    mct2 = scdb_mctnetwork("test")

    return(paste0("Recomputed network flow model is identical with sing_emb_wt10: ",
           identical(mct1@network,mct2@network)))
}