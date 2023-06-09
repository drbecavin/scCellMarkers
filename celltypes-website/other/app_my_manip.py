import json
from flask import Flask, jsonify, render_template, request, Markup
import graphviz
import io
import pandas as pd

app = Flask(__name__)

# Load graph from .dot file
graph = graphviz.Source.from_file("graph/endothelial.dot")
table = pd.read_csv("test_table.csv")    

#def print_rows_table(node_name):
#    x = table[table["celltype"] == node_name]
#    print(x)
#    return(x)


from datetime import date

# @app.route('/node_info')
@app.route('/node_info')
def node_info():
    print("node_info")
    node_name = request.args.get("node_name")

    # Look up the node info based on the node name
    #    node_info = my_graph.get_node_info(node_name)
    node_info = "Info about the node"

    # Get the rows of the table corresponding to the node name
    x = table[table["celltype"] == node_name]
    rows = x.to_dict(orient='records')
    # Return the node name and info as a JSON object
    #return json.dumps({
    #    'name': node_name,
    #    'info': node_info,
    #    'table_rows': rows
    #})
    return jsonify(x)
    
    # original line to return name:
    # return jsonify(node_name)
    # return jsonify(node_name), node_info
    
    # tried the following two (not together) but did not work
    # return jsonify(node_name, dict_info)
    # return jsonify(dict_info)
    
    # tried the following code that did not work
    # resp = jsonify(dict_info)
    # resp.status_code = 200
    # print("this is resp", resp)
    # return resp
    
    # test to see if the problem is the dictionary or the 
    # values
    # dict_test = {"x": 24356, "y": 3456}
    # resp = jsonify(dict_test)
    # resp.status_code = 200
    # print(resp)
    # return resp
    
    # test to see if the problem is the dictionary or the 
    # values
    # dict_test = {"x": 24356, "y": 3456}
    # resp = jsonify(dict_test)
    # print(resp)
    # return resp
    
    #return jsonify(
    #    name = node_name, 
    #    info = node_info
    #)
    # return dict_test

@app.route("/")
def index():
    
    graph_svg = graph.render(format="svg")

    # Convert the SVG image to a text file like object
    svg_text = ""
    with open(graph_svg, "r") as file:
        for line in file:
            svg_text+=line
    
    # Pass the file-like object to the template
    return render_template("index_my_manip.html", graph_svg=svg_text)

if __name__ == "__main__":
    app.run(debug = True)