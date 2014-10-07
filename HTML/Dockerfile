FROM ubuntu

RUN apt-get update
RUN apt-get -y install python-software-properties software-properties-common
RUN add-apt-repository ppa:j-4/vienna-rna
RUN apt-get update
RUN apt-get -y install gbrowse 
RUN apt-get -y install default-jdk 
RUN apt-get -y install vienna-rna wget

RUN mkdir -p /tmp
RUN cd /tmp && wget http://mn.eng.hawaii.edu/~garmire/weka.jar && cp weka.jar /usr/lib/jvm/default-java/jre/lib/ext/.

RUN apt-get -y install build-essential
RUN cd /tmp && wget http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz && tar zxvf miRanda-aug2010.tar.gz && cd miRanda-3.3a && ./configure && make && make install

RUN apt-get -y install r-base r-base-dev

RUN mkdir -p /tmp/test
RUN cd /tmp && wget http://mn.eng.hawaii.edu/~garmire/MirMarkRNApl.tgz && tar zxvf MirMarkRNApl.tgz && cp MirMarkRNApl/*.pl /usr/share/perl5
RUN cd /tmp && wget http://mn.eng.hawaii.edu/~garmire/MirMark.tgz && tar zxvf MirMark.tgz && cp MirMark/*.* /usr/local/bin
