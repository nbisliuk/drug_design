import sqlite3

# Create connection to a database
conn = sqlite3.connect('chembl_28/chembl_28.db')

# Create a cursor object
c = conn.cursor()

c.execute("""
    SELECT DISTINCT canonical_smiles FROM compound_structures WHERE molregno IN (
        SELECT molregno FROM activities WHERE
            standard_value IS NOT NULL AND
            standard_units = 'nM' AND
            standard_type IN ("IC50", "EC50", "Ki", "Kb", "Kd")
        INTERSECT
        SELECT molregno FROM molecule_dictionary WHERE molecule_type = 'Small molecule'
    )
""")

items = c.fetchall()

# Print fetched results
for item in items[:10]:
    print(item)

# Save the results into file
with open('data/canonical_smiles.smi', "w") as output:
    output.writelines(["%s\n" % item for item in items])

# Close the connection
conn.close()