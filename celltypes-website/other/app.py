from flask import Flask, jsonify, render_template, request, Markup
import graphviz
import io

app = Flask(__name__)

# Load graph from .dot file
graph = graphviz.Source.from_file("graph/endothelial.dot")
    

def get_nodename(node_name):
    print(node_name)


@app.route("/node_click")
def node_click():
    node_name = request.args.get("node_name")

    
    ####
    
    result = get_nodename(node_name)
    
    ####
    
    
    return jsonify(result)

@app.route('/node_info')
def node_info():
    print("node_info")
    node_name = request.args.get("node_name")

    # Look up the node info based on the node name
    #    node_info = my_graph.get_node_info(node_name)
    node_info = "Info about the node"

    dict_info = {'name': node_name, 'info': node_info}
    # Return the node name and info as a JSON object
    return jsonify(node_name)

@app.route("/")
def index():
    
    graph_svg = graph.render(format="svg")

    # Convert the SVG image to a text file like object
    svg_text = ""
    with open(graph_svg, "r") as file:
        for line in file:
            svg_text+=line
    
    # Pass the file-like object to the template
    return render_template("index.html", graph_svg=svg_text)

if __name__ == "__main__":
    app.run()