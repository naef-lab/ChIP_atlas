CREATE DATABASE peak_database;

USE peak_database;

-- Define table structure as per your text file structure: chr start   end gene    score   strand  promoter_id antigen exp_id  rel_start   rel_end
CREATE TABLE peak_table (
    chr VARCHAR(255),
    start_ INT,
    end_ INT,
    gene VARCHAR(255),
    score FLOAT,
    strand VARCHAR(255),
    promoter_id VARCHAR(255),
    antigen VARCHAR(255),
    exp_id VARCHAR(255),
    rel_start INT,
    rel_end INT,
    -- Define columns as per your text file structure
);

-- Load data from text file into table for loop in results/mm10/Peak_tables/Window_pm5kb/*.bed


LOAD DATA LOCAL INFILE '/path/to/your/file.csv'
INTO TABLE my_table
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(column1, column2, column3, column4);