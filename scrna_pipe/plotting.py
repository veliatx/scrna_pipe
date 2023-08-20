import requests


def plot_string_network(gene_list):
    """
    """
    string_api_url = "https://version-11-5.string-db.org/api"

    output_format = "highres_image"
    method = "network"

    params = {

        "identifiers" : "\r".join(gene_list),
        "species" : 9606, 
        "limit" : 1, 
        "echo_query" : 1, 
        "caller_identity" : "www.awesome_app.org" 

    }

    request_url = "/".join([string_api_url, output_format, method])

    response = requests.post(request_url, data=params)

    return response.content