#!/bin/python
#Jan 2018 - Kuan-Lin Huang @ WashU - 
# Kuan's first script in 2018!!!
 
import sys
import getopt

def main():
    def usage():
        print """
    USAGE: convert_pml_musite_color.py [-h] <clusterID> <PyMol File> 
     <filename>    input file
        """

    #use getopt to get inputs
    # try:
    #     opts, args = getopt.getopt(sys.argv[1:], 'h') #:after option meaning required arguments
    # except getopt.GetoptError:
    #     print "liftover_CharGer_result.py <CharGerSummaryFile> <inputFile>"

    # for opt, arg in opts: #store the input options
    #     if opt == '-h': # h means user needs help
    #         usage(); sys.exit()

    # file IO
    args = sys.argv[1:]
    if len(args) != 1:
        usage(); sys.exit("input file missing")

    try:
        pymolF = open(args[0],"r")
    except IOError:
        print("PyMol File , args[1], does not exist!")

    outFH = args[0].replace(".pml","") + ".colored.pml"
    outF = open(outFH, "w")
    sys.stdout = outF

    # colors
    reds = ["red","tv_red","raspberry","darksalmon","salmon","deepsalmon","warmpink","firebrick","ruby","chocolate","brown"]
    blues = ["blue","tv_blue","marine","slate","lightblue","skyblue","purpleblue","deepblue","density"]
    magentas = ["magenta","lightmagenta","hotpink","pink","lightpink","dirtyviolet","violet","violetpurple","purple","deeppurple"]
    ptm_i = 0
    mut_i = 0

    #read pymol file
    for line in pymolF:
        line=line.strip()

        # lines to print the residues
        if line.startswith("sele"):
            # only print variant if its in the variant class as we recorded in the cluster file
            F = line.split(" ")
            color_var = F[1]
            color_var_split = color_var.split("_")
            #originalColor = color_var_split[0]
            originalColor = F[8].replace(",","")
            var = color_var_split[2]
            
            varColor = "forest" # default to green, check for error when this happens
            # print color_var
            if "ptm" in var:
                varColor = reds[ptm_i]
                ptm_i = ptm_i + 1
            else:
                varColor = blues[mut_i]
                mut_i = mut_i + 1

            if ptm_i == len(reds):
                ptm_i = 0
            if mut_i == len(blues):
                mut_i = 0

            updatedLine = line.replace(originalColor, varColor) # doesn't change grey of the centroid for now
            updatedLine = updatedLine.replace("grey", varColor)
            print updatedLine

        # lines to adjust the backbone chain
        elif line.startswith("color"):
            F = line.split(" ")
            F[1] = "grey90,"
            print " ".join(F)
	
	# last line
	elif line.startswith("#mplay"):
	    print line
	    print "#png " + outFH + ".png ,width=1200,height=1200, dpi=300, ray=1;"
	    print "#quit"

        # set-ups
        else:
            print line
            # print an extra line for background
            if line.startswith("bg_color"):
                print "set ray_opaque_background, off;" # change background color
        

    pymolF.close()



if __name__ == "__main__":
    main()
