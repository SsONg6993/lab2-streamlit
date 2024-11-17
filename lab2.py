import streamlit as st
import requests
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "52a396ee4431f8ed1f73b04cb5d9373d",
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,
        "searchbiogridids": True,
        "includeInteractors": True
    }
    response = requests.get(biogrid_url, params=params)
    network = response.json()
    network_df = pd.DataFrame.from_dict(network, orient='index')
    network_df.OFFICIAL_SYMBOL_A = [gene.upper() for gene in network_df.OFFICIAL_SYMBOL_A]
    network_df.OFFICIAL_SYMBOL_B = [gene.upper() for gene in network_df.OFFICIAL_SYMBOL_B]
    return network_df


def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606
    }
    response = requests.get(string_url, params=params)
    network = response.json()
    network_df = pd.json_normalize(network)
    network_df['preferredName_A'] = network_df['preferredName_A'].str.upper()
    network_df['preferredName_B'] = network_df['preferredName_B'].str.upper()
    return network_df


def generate_network(dataframe):
    network_graph = nx.Graph()
    if 'OFFICIAL_SYMBOL_A' in dataframe.columns and 'OFFICIAL_SYMBOL_B' in dataframe.columns:
        network_graph = nx.from_pandas_edgelist(dataframe, "OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B")
    elif 'preferredName_A' in dataframe.columns and 'preferredName_B' in dataframe.columns:
        network_graph = nx.from_pandas_edgelist(dataframe, "preferredName_A", "preferredName_B")
    else:
        st.warning("PPI data not found!")
    return network_graph


def get_centralities(network_graph):
    degree_centrality = nx.degree_centrality(network_graph)
    betweenness_centrality = nx.betweenness_centrality(network_graph)
    closeness_centrality = nx.closeness_centrality(network_graph)
    eigenvector_centrality = nx.eigenvector_centrality(network_graph, 200)
    pagerank_centrality = nx.pagerank(network_graph)

    return {
        'Degree Centrality': degree_centrality,
        'Betweenness Centrality': betweenness_centrality,
        'Closeness Centrality': closeness_centrality,
        'Eigenvector Centrality': eigenvector_centrality,
        'PageRank Centrality': pagerank_centrality
    }

def interpret_shared_top_nodes(centralities):
    # Retrieve the top 5 nodes for each centrality measure
    top_5_nodes_by_measure = {}
    for measure, values in centralities.items():
        top_5_nodes = sorted(values.items(), key=lambda x: -x[1])[:5]
        top_5_nodes_by_measure[measure] = {node for node, _ in top_5_nodes}
    
    # Find common nodes across all measures
    common_top_nodes = set.intersection(*top_5_nodes_by_measure.values())
    
    if common_top_nodes:
        st.write("### Common Top Nodes Across All Measures")
        st.write(f"These nodes are ranked in the top 5 across all measures: {', '.join(common_top_nodes)}")
        st.write("These nodes likely play critical roles in the network, showing high interaction and connectivity.")
    else:
        st.write("### No Common Top Nodes Across All Measures")
        st.write("There are no nodes that appear in the top 5 across all centrality measures.")

    return common_top_nodes


st.title('LAB 2 - ONG SHUN SHENG')
st.markdown("### Enter Protein Information Below:")

protein_id = st.text_input("Enter Protein ID:")
database_choice = st.selectbox("Choose Database", ["BioGRID", "STRING"])
retrieve = st.button('Retrieve')

if retrieve:
    if database_choice == "BioGRID":
        ppi_data = retrieve_ppi_biogrid(protein_id)
    else:
        ppi_data = retrieve_ppi_string(protein_id)

    if not ppi_data.empty:
        col1, col2 = st.columns(2)

        network_graph = generate_network(ppi_data)
        slayout = nx.spring_layout(network_graph, seed=123)

        node_size = 20 if database_choice == "BioGRID" else 1000

        with col1:
            st.subheader("PPI Data Information")
            st.dataframe(ppi_data)
            st.write(f"**Number of edges:** {network_graph.number_of_edges()}")
            st.write(f"**Number of nodes:** {network_graph.number_of_nodes()}")

            st.markdown("### Network Visualization")
            plt.figure(figsize=(8, 6)) 
            nx.draw(network_graph, slayout, with_labels=False, node_size=node_size, node_color='lightblue')
            st.pyplot(plt)

        with col2:
            st.subheader("Centrality Measures")
            centralities = get_centralities(network_graph)

            for measure, values in centralities.items():
                st.write(f"**{measure}:**")
                top_5_nodes = sorted(values.items(), key=lambda x: -x[1])[:5]
                for index, (node, value) in enumerate(top_5_nodes, start=1):
                    st.write(f"{index}. {node}: {value:.4f}")

                top_5_node_names = [node for node, _ in top_5_nodes]

                if measure == 'Degree Centrality':
                    st.write(f"   - Proteins '{top_5_node_names}' have higher Degree Centrality, indicating they will interact with many other proteins, suggesting they may be a hub proteins. Hub proteins are often crucial in cellular processes and are more likely to be essential for survival or associated with disease pathways.")

                elif measure == 'Betweenness Centrality':
                    st.write(f"   - Proteins '{top_5_node_names}' have higher Betweenness Centrality, indicating they often serve as connectors or regulators of information flow between different parts of the network, which can make them important for maintaining network integrity and facilitating communication between modules or clusters.")

                elif measure == 'Closeness Centrality':
                    st.write(f"   - Proteins '{top_5_node_names}' have higher Closeness Centrality, indicating they influence many parts of the network efficiently and are often involved in signal transduction or regulatory functions")

                elif measure == 'Eigenvector Centrality':
                    st.write(f"   - Proteins '{top_5_node_names}' have higher Eigenvector Centrality, indicating they are the influential proteins that are connected to other influential proteins, potentially marking out proteins that play central roles within highly active or essential functional regions of the network.")

                elif measure == 'PageRank Centrality':
                    st.write(f"   - Proteins '{top_5_node_names}' have PageRank Centrality, indicating they are well-connected and linked to other prominent proteins, which helping to highlight functionally essential or central proteins in signaling pathways or disease-related clusters.")

                plt.figure(figsize=(8, 6))
                nx.draw(network_graph, slayout, with_labels=False, node_size=node_size, node_color='lightblue')
                nx.draw_networkx_nodes(network_graph, slayout, nodelist=top_5_node_names, node_size=node_size*5, node_color='orange')  
                st.pyplot(plt)


        common_top_nodes = interpret_shared_top_nodes(centralities)
        
        st.success("Data visualized successfully!")
    else:
        st.error("No PPI data found for the provided protein ID.")
