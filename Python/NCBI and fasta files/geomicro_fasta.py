from Bio import SeqIO
import re


def parse_genbank(gb_file):
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

            
        acc_num = gb_annotations["db_source"]
        acc_nm = extract(acc_num, r'\w+\.\w+')
        location_by_source = index_genbank_features(gb_record,"source","country")


        print acc_nm+" "+location_by_source
             
             
             





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
