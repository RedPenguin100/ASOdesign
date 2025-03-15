__All files may also be found in the drive, under Project/Model/Human genome v34__

In this file, all files(except DBs) we expect to find were downloaded from https://www.gencodegenes.org/human/release_34.html and are found in the google drive.

For the DBs, only the google-drive / manual creation can help. Running the `read_human_genome.py` will automatically create the DB, but it is very slow.

### First step
Download all files from the drive under "/Project/Model/Human genome v34" and put them in this directory
### Second step
Extract all files from the drive, using(in Linux) 
`gunzip -k <file_path>`
In Windows you may try to use gzip from command line.
### Test
Run the `read_human_transcriptome.py` for a few second, see if it starts up.

Enjoy!
