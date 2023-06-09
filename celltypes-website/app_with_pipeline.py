from flask import Flask, jsonify, render_template, request
import graphviz
import load
import markers_research 


app = Flask(__name__)

# Load graph from .dot file
graph = graphviz.Source.from_file(
    load.paths()['DATA_PATH'] + "/graphviz_tree/trees_V8/endothelial/Digraph.gv.dot")
print(load.paths()['DATA_PATH'])


# Define function to change celltype compartment
# e.g. immune or stromal
def choose_dot_file(compart):
    print(load.paths()[
        'DATA_PATH'] + f"/graphviz_tree/trees_V8/{compart}/Digraph.gv.dot")
    return graphviz.Source.from_file(load.paths()[
        'DATA_PATH'] + f"/graphviz_tree/trees_V8/{compart}/Digraph.gv.dot")


@app.route('/node_info')
def node_info():
    print("node_info")
    node_name = request.args.get("node_name")
    return jsonify(node_name)


@app.route('/get_filtered_table')
def get_filtered_table():
    print("get_filtered_table")
    node_name = request.args.get("node_name")

    # Markers research
    filtered_df = markers_research.search_cell(node_name, n=50)

    tempfile_name = f"static/tmp/{node_name}.tsv"

    # Convert the filtered dataframe to HTML
    filtered_df_html = f"<a href=\"{tempfile_name}\">Download_table</a><br><br>"

    # Convert the filtered dataframe to HTML
    filtered_df_html += filtered_df.to_html(index=False)

    return jsonify(filtered_df_html)


@app.route("/", methods=["GET", "POST"])
def index():

    selected_compart = "endothelial"

    if request.method == "POST":
        selected_compart = request.form["compart"]
        # use the selected option here
        print(f"The selected compartment is {selected_compart}")
        global graph
        graph = choose_dot_file(selected_compart)

    graph_svg = graph.render(format="svg")

    # Convert the SVG image to a text file like object
    svg_text = ""
    with open(graph_svg, "r") as file:
        for line in file:
            svg_text += line.replace('<svg ', '<svg id="graph" ')
    
    # Pass the file-like object to the template
    return render_template("index_with_pipeline.html", graph_svg=svg_text,
                           selected_compart=selected_compart)


if __name__ == "__main__":
    app.run()
