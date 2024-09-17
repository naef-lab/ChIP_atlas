import os
import pandas as pd
from sqlalchemy import create_engine

# Set up MySQL connection
db_user = 'jbreda'
db_password = 'password'
db_host = 'localhost'
db_port = '3306'
db_name = 'Chip_Atlas_Peaks_DB'

# Create the MySQL engine
engine = create_engine(f"mysql+mysqlconnector://{db_user}:{db_password}@{db_host}/{db_name}")

# Path where your text files are located
Genome = ['hg38','mm10']
for genome in Genome:
    table_name = f'Peak_table_{genome}'
    file_path = f'results/{genome}/Peak_tables/Window_pm5kb/'
    # Loop through all the text files in the directory
    for filename in os.listdir(file_path):
        if filename.endswith('.bed'):  # or other text file extensions
            file_full_path = os.path.join(file_path, filename)
            
            # Read the CSV file into a Pandas DataFrame
            df = pd.read_csv(file_full_path,sep='\t')
            
            # Write the DataFrame to the MySQL table
            df.to_sql(name=table_name, con=engine, if_exists='append', index=False)
            print(f"File {filename} imported successfully.")
    print(f"All files for {genome} imported successfully.")
print("All files imported successfully.")
