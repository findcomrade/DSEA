# install.packages("RMySQL")
library(RMySQL)

# Connection:
mydb = dbConnect(MySQL(), user='dbychkov', password='', dbname='chembl_16', host='lindex.giu.fi')

# Listing Tables and Fields:
dbListTables(mydb)
dbListFields(mydb, 'compound_records')

# Retrieving data:
rs = dbSendQuery(mydb, "SELECT molregno FROM atc_classification WHERE who_name = 'Dasatinib'")
data = fetch(rs, n=-1)