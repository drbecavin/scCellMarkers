from flask import Flask, jsonify, render_template, request, Markup
import graphviz
import io
import pandas as pd

app = Flask(__name__)

# Load graph from .dot file
graph = graphviz.Source.from_file("graph/endothelial.dot")
    
# Define your pandas dataframe here
df = pd.DataFrame({
        'celltype': ['Endothelial', 'Blood vessels EC', 'Lymphatic EC', 'EC arterial', 'EC capillary',
                     'EC venous', 'Lymphatic EC differentiating', 'Lymphatic EC mature',
                     'EC aerocyte capillary', 'EC general capillary', 'EC venous pulmonary', 'EC venous systemic'],
        'x': [25, 30, 35, 40, 45, 25, 30, 35, 40, 45, 11, 12],
        'y': ['M', 'M', 'F', 'M', 'F', 'M', 'M', 'F', 'M', 'F', 'M', 'F']
    })

def get_nodename(node_name):
    print(node_name)


@app.route('/node_info')
def node_info():
    print("node_info")
    node_name = request.args.get("node_name")

    return jsonify(node_name)



@app.route('/get_filtered_table')
def get_filtered_table():
    print("get_filtered_table")
    node_name = request.args.get("node_name")

    # Select rows based on the node name
    filtered_df = df[df['celltype'] == node_name]

    # Convert the filtered dataframe to HTML
    filtered_df_html = filtered_df.to_html(index=False)

    return jsonify(filtered_df_html)


@app.route("/")
def index():
    
    graph_svg = graph.render(format="svg")

    # Convert the SVG image to a text file like object
    svg_text = ""
    with open(graph_svg, "r") as file:
        for line in file:
            svg_text+=line
    
    # Pass the file-like object to the template
    return render_template("index_basic_table.html", graph_svg=svg_text)

if __name__ == "__main__":
    app.run()