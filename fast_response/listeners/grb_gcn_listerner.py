''' Script to automatically receive GCN notices for IceCube
    alert events and run followup accordingly

    Author: Alex Pizzuto
    Date:   July 2020
'''

# Download the SQLite file from the GRBweb 2 webpage
os.system("wget https://icecube.wisc.edu/~grbweb_public/GRBweb2.sqlite")

# Load the database with the sqlite3 module
db = sqlite3.connect('GRBweb2.sqlite')

Summary_table = pandas.read_sql_query("SELECT * from Summary", db)

Summary_table[Summary_table.mjd==max(Summary_table.mjd)]