import sqlite3 as db

class GeneQuery2:
    def __init__(self):
        self.connection = db.connect('database/GeneQuery2.db')

    def get_clusters(self):
        cursor = self.connection.cursor()
        find_clusters = cursor.execute('SELECT name from sqlite_master where type="table"')
        clusters_list = [n[0] for n in cursor.fetchall()]
        return clusters_list



