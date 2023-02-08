import pwd
import os
import base64
from io import BytesIO
#from web_app.functions import plot_umap_gene

import pandas as pd

from flask import Flask, render_template, request, flash, redirect, url_for, send_from_directory, session, send_file, make_response, render_template_string
from werkzeug.utils import secure_filename

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from functions import read_h5ad, plot_umap, plot_umap_gene_violin, plot_dotplot, analysis_clustering, diff_exp_table

UPLOAD_FOLDER = 'uploads/'
ALLOWED_EXTENSIONS = {'csv', 'tsv', 'mtx', 'h5ad'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'secret'

adata=None
adata2=None

## Homepage
@app.route('/')
def homepage():
    return render_template('homepage.html') 


## Upload
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            print(type(file.filename))
            session['file_name']=file.filename
            flash('File downloaded sucsessfully')
            print('File downloaded')
            return redirect(url_for('download_file', name=filename))
    return render_template("upload.html")
		
@app.route('/uploads/<name>')
def download_file(name):
    #flash('File downloaded sucsessfully')
    return send_from_directory(app.config["UPLOAD_FOLDER"], name)

@app.route("/uploads/forward/", methods=['POST'])
def move_forward():
    #Moving forward code
    forward_message = "Reading..."
    global adata
    adata=read_h5ad(f"uploads/{session['file_name']}")
    print_adata=adata.__repr__
    return render_template('upload.html', forward_message=forward_message, print_adata=print_adata);


## Analysis Page
@app.route("/analysis")
def analysis():
    print(session['file_name'])
    global adata
    adata=read_h5ad(f"uploads/{session['file_name']}")
    print_adata=adata.__repr__
    global adata2
    adata2=read_h5ad('./uploads/clustering_results.h5ad')
    print_adata2=adata2.__repr__
    return render_template("analysis.html", print_adata=print_adata,
                            print_adata2=print_adata2)

# Analyze data
@app.route("/analysis/forward/")
def analysis_forward():
    if request.args:
        n_dict = dict(request.args)
        n = pd.to_numeric(n_dict.get('n_neighbors'))
        if n:
            print(n, n_dict)
        else:
            n = 10
            print('No argument', n, n_dict)
    else:
        n = 10
        print('No input')
    global adata
    print_adata=adata.__repr__
    global adata2
    #adata2 = adata.copy()
    adata2=analysis_clustering(adata, n_neighbors=n)
    forward_message = "Analysis done"
    print_adata2=adata2.__repr__
    adata2.write('./uploads/clustering_results.h5ad')
    return render_template('analysis.html', print_adata=print_adata, 
                            forward_message=forward_message, print_adata2=print_adata2)

# Save analysis results
@app.route("/analysis/save", methods=['GET'])
def save_file():
    global adata
    print_adata=adata.__repr__
    global adata2
    print_adata2=adata2.__repr__
    try: 
        adata2.write('./uploads/clustering_results.h5ad')
        print('Saving analysis file')
        return send_file('./uploads/clustering_results.h5ad')
    except:
        forward_message = "Could not save the file"
        print('Could not save the file')
        return render_template('analysis.html', print_adata=print_adata, 
                            forward_message=forward_message, print_adata2=print_adata2)


@app.route("/results/direct", methods=['GET'])
def directly_to_results():
    global adata
    print_adata=adata.__repr__
    global adata2
    adata.write('./uploads/clustering_results.h5ad')
    return render_template('analysis.html', print_adata=print_adata 
                            )


# Results page
@app.route("/results")
def results():
    global adata2
    adata2=read_h5ad('./uploads/clustering_results.h5ad')

    maps_list = list(adata2.obsm)[1:]
    gene_name_list = list(adata2.var_names)
    print(maps_list)

    # Get input from user 
    if request.args:
        dict_arg = dict(request.args)
        print(dict_arg)

        if dict_arg.get('map'):
            global map
            map=dict_arg.get('map')
            print (map)
        elif dict_arg.get('gene'):
            global gene
            gene=dict_arg.get('gene')
            print (gene)
        else:
        # form_dict = dict(request.args)
        # form = form_dict.get('map')
            print('Something is wrong with request.args')

    else:
        map =maps_list[0]
        gene =gene_name_list[0]
        print('No requested args, default setting')
        print(map, gene)

    try:
        #/ Plot clusters
        fig1 = plot_umap(adata2, map=map)
        fig1.savefig('./uploads/clustering_results.pdf')
        buf = BytesIO()
        fig1.savefig(buf, format="png") # Save it to a temporary buffer.
        img_data_l = base64.b64encode(buf.getbuffer()).decode("ascii") # Embed the result in the html output.
        #return f"<img src='data:image/png;base64,{data}'/>"
        plt.close()
        buf.close()

        #/ Plot genes
        fig = plot_umap_gene_violin(adata2, map=map, color=gene)
        fig.savefig('./uploads/gene_umap.pdf')
        buf = BytesIO()
        fig.savefig(buf, format="png")
        img_data_g = base64.b64encode(buf.getbuffer()).decode("ascii")
        plt.close()
        buf.close()

        # Make table 
        table=diff_exp_table(adata2)
        table.to_csv('./uploads/DiffGene_table.csv')

        # Plot dotplot
        fig2 = plot_dotplot(adata2)
        fig2.savefig('./uploads/DiffGene_dotplot.pdf')
        buf = BytesIO()
        fig2.savefig(buf, format="png")
        img_data_d = base64.b64encode(buf.getbuffer()).decode("ascii")
        plt.close()
        buf.close()
        return render_template("results.html", 
                            data1=img_data_l, maps=maps_list,
                            data2=img_data_g, gene_names=gene_name_list,
                            data3=img_data_d
                          )
    except:
        print('Something went wrong')
        return render_template("results.html"
                          )

# Save figures
@app.route("/results/save_umap", methods=['GET'])
def save_umap():
    try: 
        print('Saving umap')
        return send_file('./uploads/clustering_results.pdf')
    except:
        forward_message = "Could not save the file"
        print('Could not save the file')
        return render_template("results.html")

@app.route("/results/save_dotplot", methods=['GET'])
def save_dotplot():
    try: 
        print('Saving dotplot')
        return send_file('./uploads/DiffGene_dotplot.pdf')
    except:
        forward_message = "Could not save the file"
        print('Could not save the file')
        return render_template("results.html")

@app.route("/results/save_gene_umap", methods=['GET'])
def save_gene_umap():
    try: 
        print('Saving umap')
        return send_file('./uploads/gene_umap.pdf')
    except:
        forward_message = "Could not save the file"
        print('Could not save the file')
        return render_template("results.html")

@app.route("/results/save_table", methods=['GET'])
def save_table():
    try: 
        print('Saving table')
        return send_file('./uploads/DiffGene_table.csv')
    except:
        forward_message = "Could not save the file"
        print('Could not save the file')
        return render_template("results.html")


# Show tables
@app.route("/tables")
def show_tables():
    # global adata2
    # adata2=read_h5ad('./uploads/clustering_results.h5ad')
    #try:
    # table=diff_exp_table(adata2)
    # print(table)
    # df_html = table.to_html()
    #resp = make_response(render_template_string(df_html))
    #resp = make_response(render_template('tables.html',table=table))
    #return resp

    return render_template("tables.html")
    # except:
    #     print('Something went wrong')
    #     return render_template("tables.html"
    #                       )




if __name__ == '__main__':
    app.run(debug=True)
