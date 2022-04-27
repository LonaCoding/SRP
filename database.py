import sqlite3 as db


class GeneQuery:
    def __init__(self, database_number: int):
        self.database_number = database_number
        self.connection = db.connect(f'/local/www/htdocs/webapp/database/GeneQuery{self.database_number}.db')
    def get_clusters(self) -> list:
        cursor = self.connection.cursor()
        find_clusters = cursor.execute('SELECT name from sqlite_master where type="table"')
        clusters_list = [n[0] for n in cursor.fetchall()]
        return clusters_list

    def find_genes(self, cluster: str, gene: str) -> list:
        cursor = self.connection.cursor()

        find = cursor.execute(f'SELECT * from {cluster} WHERE gene = "{gene}"')
        result = cursor.fetchall()
        return result

    def gene_counts(self, cluster: str):
        cursor = self.connection.cursor()
        find = cursor.execute(f'SELECT COUNT(*) FROM {cluster}')
        result = cursor.fetchall()
        return result


class FoundGene2:
    def __init__(self, row: list):
        self.name = row[0][0]
        self.p_value = row[0][-1]
        self.log2FC = row[0][2]

class FoundGene1:
    def __init__(self, row: list):
        self.name = row[0][0]
        self.p_value = row[0][-1]
        self.log2FC = row[0][2]



class GeneCounts(GeneQuery):
    def __init__(self,database_number):
        super().__init__(database_number)
        self.genes_per_cluster = [self.gene_counts(n)[0][0] for n in self.get_clusters()]
        self.total_genes = sum(self.genes_per_cluster)
