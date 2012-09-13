from Bio import SeqIO
import re

#Sets up a fasta to help create a biogeography tree
#Use this on a genbank file to extract the acceession number, location,
#and isolarion (basically what the sample was taken from)
#Comments coming soonish
#


def parse_genbank(gb_file):

    out_geo = open("new"+"_"+gb_file+".fasta","w+")
    out_geo_ref = open("ref"+"_"+gb_file+".csv","w+")
    
    recs = SeqIO.parse(gb_file,"genbank")
    for gb_record in recs:
        
        try:
            gb_feature = gb_record.features[2]
        except IndexError:
                gb_feature = gb_record.features[1]

        def index_genbank_features(gb_record, feature_type, qualifier) :
            answer = dict()
            for (index, feature) in enumerate(gb_record.features) :
                if feature.type==feature_type :
                    if qualifier in feature.qualifiers :
                        for value in feature.qualifiers[qualifier] :
                            if value in answer :
                                print ("WARNING - Duplicate key %s for %s \
                                          features %i and %i") \
                                          % (value, feature_type, answer[value],
                                          index)
                            else :
                                answer[value] = index
                return answer
        location_by_source = index_genbank_features(gb_record,"source","country")
        isolation_by_source = index_genbank_features(gb_record,"source","isolation_source")

        location = extract(location_by_source,r'\s[^0-9:"]+')
        source = extract(isolation_by_source,r'[\w+\s]+')

        #print gb_record.id
        #print location
        #print source
        #print gb_record.seq

        loc = str(location)

        out_geo.write(">"+gb_record.id+" "+location+"\n"+str(gb_record.seq)+"\n")
        out_geo_ref.write(gb_record.id+","+location+","+source+"\n")

        
    out_geo.close()
    out_geo_ref.close()

#####
#
#This function sets up the reguar expression routine for pulling out information
#in the fasta file. This will convert the extracted pattern to a string
#
#####
def extract(src, pattern):
     src_str = str(src)
     matchObj = re.search(pattern, src_str)
     if matchObj:
          return matchObj.group()
     else:
          return "NONE"

#Replace the stuff between the quotes with your file name
parse_genbank("Northup_16S_lava.gb")             
             
             






