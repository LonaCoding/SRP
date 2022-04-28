from flask import Flask
from flask import render_template
from flask import request
import sqlite3
import os
import datetime as dt
import database
import template_inserts

app = Flask(__name__)

date = dt.date.today().year  # Date for the footer

''' App routes and functions for relevant webpages'''


@app.route('/')  # Function and route for the homepage
def index():
    with open('/local/www/htdocs/webapp/static/text/main_page/introduction.txt', 'r', encoding='utf-8') as file:
        introduction = file.read()
    with open('/local/www/htdocs/webapp/static/text/main_page/project_overview1.txt', 'r', encoding='utf-8') as file:
        project_overview1 = file.read()
    with open('/local/www/htdocs/webapp/static/text/main_page/project_overview2.txt', 'r', encoding='utf-8') as file:
        project_overview2 = file.read()
    return render_template('index.html', introduction=introduction, project_overview1=project_overview1,
                           project_overview2=project_overview2, date=date)


@app.route('/pipeline/<int:number>/')  # Function and route for pipeline 1 and 2. <int> changes depending on the link
# the user clicks
def pipeline(number):
    pipeline_num = number

    # Process relevant text and image files for pipeline page
    template_insert = template_inserts.TemplateInsert(number)
    texts = template_insert.get_text()

    figure_legends = texts[0]
    figure_paragraphs = texts[1]
    methods_results = texts[2]

    figures = template_insert.get_images()
    return render_template('analysis_pipelines.html', pipeline_num=pipeline_num, figure_legends=figure_legends,
                           figure_paragraphs=figure_paragraphs, figures=figures, methods_results=methods_results,
                           zip=zip, type=type, str=str, len=len,
                           list=list, date=date)  # All relevant variables required to be parsed into the html page
    # via Jinja2


@app.route('/pipeline/<int:number>/query', methods=['GET', 'POST'])  # Route and function for the gene query webpage
# where the user can search the database for genes from a given cluster. This is done via POST request
def gene_query(number):
    pipeline_num = number

    selected_database = database.GeneQuery(pipeline_num)
    clusters = sorted(selected_database.get_clusters())
    gene_count_info = database.GeneCounts(pipeline_num)
    total_genes = gene_count_info.total_genes

    return render_template('genequery.html', pipeline_num=pipeline_num, clusters=clusters, total_genes=total_genes,
                           len=len, date=date)


@app.route('/pipeline/<int:number>/query/search', methods=['POST'])  # Route and function that recieves the POST
# requests and executes a SQL command to search a given database
def search_database(number):
    pipeline_num = number
    gene = request.form.get('gene').upper()
    cluster = request.form.get('cluster')

    selected_database = database.GeneQuery(pipeline_num)
    find_gene = selected_database.find_genes(cluster, gene)

    if len(find_gene) == 0:
        found = False
        results = f'{gene} could not be found in {cluster}'
    else:
        found = True
        if pipeline_num == 1:
            results = database.FoundGene1(find_gene)
        else:
            results = database.FoundGene2(find_gene)

    return render_template('genesearch.html', pipeline_num=pipeline_num, gene=gene, cluster=cluster, found=found,
                           results=results, str=str(), date=date)


@app.route('/references')  # References page
def references():
    return render_template('References.html', date=date)


if __name__ == '__main__':
    app.run(debug=True)
