MirMark: a site-level and UTR-level classifier for miRNA target prediction
==========================================================================
   
Installation
------------

### Linux

First, if you haven't already, install all the dependencies:

    sudo add-apt-repository ppa:j-4/vienna-rna
    sudo apt-get update

    sudo apt-get -y install python-software-properties software-properties-common
    sudo apt-get -y install gbrowse 
    sudo apt-get -y install default-jdk 
    sudo apt-get -y install vienna-rna wget

    cd /tmp && wget http://mn.eng.hawaii.edu/~garmire/weka.jar && sudo cp weka.jar /usr/lib/jvm/default-java/jre/lib/ext/.

    sudo apt-get -y install build-essential

    cd /tmp && wget http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz && tar zxvf miRanda-aug2010.tar.gz && cd miRanda-3.3a && ./configure && make && sudo make install

    cd /tmp && wget http://mn.eng.hawaii.edu/~garmire/MirMarkRNApl.tgz && tar zxvf MirMarkRNApl.tgz && sudo cp MirMarkRNApl/*.pl /usr/share/perl5
    cd /tmp && wget http://mn.eng.hawaii.edu/~garmire/MirMark.tgz && tar zxvf MirMark.tgz && sudo cp MirMark/*.pl /usr/local/bin

Next, download this repository to local machine:

    git clone --depth 1 https://github.com/lanagarmire/MirMark.git

After downloading the repository, enter the `Core/` folder you will see all the scripts:

    cd MirMark/Core
    ls

In order for them to run properly, copy them into your local `bin` directory:

    cp * /usr/local/bin/

Now you are good to go. You can test the script on the test files under `MirMark/Test` folder.

There are two separate scripts for predicting targets in site and UTR level, respectively.

The script for predicting targets in site level

    siteFeaturesARFF.pl <mir-fasta-file> <utr-fasta-file> <phastcon-utr-scores> <pair-file> <output-prefix> MirMark/Test/rf.site.model

For `<mir-fasta-file>`, you can retrieve it from mirbase.org. For `<utr-fasta-file>`, go to the UCSC table browser.

The format for `<phastcon-utr-scores>` is:

    ... see MirMark/Test/fast.txt for example

The `<Pair file>` is a TSV of miR and UTR ids, corresponding to fasta file IDs.

The `<output-prefix>` is the the output files less the extensions.

The script for predicting targets in site level

    utrFeaturesARFF.pl <mir-fasta-file> <utr-fasta-file> <phastcon-utr-scores> <pair-file> <output-prefix> MirMark/Test/rf.site.model
