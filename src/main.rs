use clap::Parser;
use netrust::{
    network::{Network, NetworkModel},
    utils::{
        load_parameters, BarabasiAlbertPars, CompletePars, ConfigurationPars, ErdosRenyiPars,
        LatticePars, NetworkPars, OutputModel, RegularPars, WattsStrogatzPars,
    },
};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub struct Args {
    #[clap(long, value_parser, default_value_t = 4)]
    pub attachment_rate: usize,
    #[clap(long, value_parser, default_value_t = 10)]
    pub average_degree: usize,
    #[clap(long, value_parser, default_value_t = 100)]
    pub degree_maximum: usize,
    #[clap(long, value_parser, default_value_t = 2)]
    pub degree_minimum: usize,
    #[clap(long, value_parser, default_value_t = false)]
    pub flag_config: bool,
    #[clap(long, value_parser, default_value_t = false)]
    pub flag_connected: bool,
    #[clap(long, value_parser, default_value = "lattice-pbc")]
    pub model_network: NetworkModel,
    #[clap(long, value_parser, default_value = "adjacency-list")]
    pub model_output: OutputModel,
    #[clap(long, value_parser, default_value_t = 100)]
    pub nxcells: usize,
    #[clap(long, value_parser, default_value_t = 100)]
    pub nycells: usize,
    #[clap(long, value_parser, default_value = "/data/")]
    pub path_target: String,
    #[clap(long, value_parser, default_value_t = 3.0)]
    pub power_law_exponent: f64,
    #[clap(long, value_parser, default_value_t = 0.0)]
    pub probability_connection: f64,
    #[clap(long, value_parser, default_value_t = 0.0)]
    pub probability_rewiring: f64,
    #[clap(long, value_parser, default_value_t = 10.0)]
    pub radius: f64,
    #[clap(long, value_parser, default_value_t = 10000)]
    pub size: usize,
    #[clap(long, value_parser, default_value_t = 5)]
    pub size_initial_cluster: usize,
}

fn main() {
    let args = Args::parse();

    let pars_net = if args.flag_config {
        match load_parameters(args.model_network) {
            Ok(params) => params,
            Err(e) => {
                eprintln!("Error loading parameters: {}", e);
                return;
            }
        }
    } else {
        match args.model_network {
            NetworkModel::BarabasiAlbert => NetworkPars::BarabasiAlbert(BarabasiAlbertPars {
                attachment_rate: Some(args.attachment_rate),
                average_degree: Some(args.average_degree),
                size: args.size,
                size_initial_cluster: Some(args.size_initial_cluster),
            }),
            NetworkModel::Complete => NetworkPars::Complete(CompletePars {
                size: args.nxcells * args.size,
            }),
            NetworkModel::Configuration => NetworkPars::Configuration(ConfigurationPars {
                degree_minimum: args.degree_minimum,
                degree_maximum: args.degree_maximum,
                power_law_exponent: args.power_law_exponent,
                size: args.size,
            }),
            NetworkModel::ConfigurationCorrelated => {
                NetworkPars::ConfigurationCorrelated(ConfigurationPars {
                    degree_minimum: args.degree_minimum,
                    degree_maximum: args.degree_maximum,
                    power_law_exponent: args.power_law_exponent,
                    size: args.size,
                })
            }
            NetworkModel::ConfigurationUncorrelated => {
                NetworkPars::ConfigurationUncorrelated(ConfigurationPars {
                    degree_minimum: args.degree_minimum,
                    degree_maximum: args.degree_maximum,
                    power_law_exponent: args.power_law_exponent,
                    size: args.size,
                })
            }
            NetworkModel::ErdosRenyi => NetworkPars::ErdosRenyi(ErdosRenyiPars {
                average_degree: Some(args.average_degree),
                probability_connection: Some(args.probability_connection),
                size: args.size,
            }),
            NetworkModel::Lattice => NetworkPars::Lattice(LatticePars {
                nxcells: args.nxcells,
                nycells: args.nycells,
            }),
            NetworkModel::LatticePBC => NetworkPars::LatticePBC(LatticePars {
                nxcells: args.nxcells,
                nycells: args.nycells,
            }),
            NetworkModel::Regular => NetworkPars::Regular(RegularPars {
                average_degree: args.average_degree,
                size: args.size,
            }),
            NetworkModel::WattsStrogatz => NetworkPars::WattsStrogatz(WattsStrogatzPars {
                average_degree: args.average_degree,
                probability_rewiring: args.probability_rewiring,
                size: args.size,
            }),
        }
    };

    let mut graph: Network;
    if args.flag_connected {
        loop {
            graph = Network::generate_network(args.model_network, pars_net.clone());
            if graph.is_connected_dfs() {
                break;
            }
        }
    } else {
        graph = Network::generate_network(args.model_network, pars_net.clone());
    }

    graph.save_to_file(&pars_net, args.model_output);
}
