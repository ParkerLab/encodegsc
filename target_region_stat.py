import base_types
import numpy
import sys
output = sys.stdout

def compute_correlation_stat(
        features,
        points_process, 
        win_size,
        num_sample = 10,
        con_percent = .99,
        block= False,
        block_coverage = .05,
        plot_type = 'line',
        aggregate_type = None,
        num_groups = 0
    ):
    """
    THis function applies the block bootstrap method in order to determine 
    the significance of a point process relative to regions 
    on a second genome.
    
    the arguments to run this are:
    -Parsed file for feature intervals
    -Parsed file for point process
    -The window size
    -The percentage for the confidence band
    -Whether or not to use the block_sampling method for the null(True to use )
    -Size of the block, if block was selected
    -Plot type: line for line plot, bar for bart plot and box for boxplot
    -aggregate_type: smoothing, either moving average or bin
    -num_group for the smoothing
    
    This function plots the desired plot 
    """
    
    chrom_index = 0
    #Iterates through each chromosome in the regions object
    for chrom_index in range(len(features.values())):
        #declares chrom_features equal to the features present in the specific chromosome
        chrom_features = features.values()[chrom_index]
        #Computes the score of the point process B in relation of the features of A
        score_array = aggregate_scores(chrom_features,points_process,win_size,chrom_features.length)
        
        #Checks if should perform a smoothing based on whether aggregate_type has a value
                    
        

        #Finds the length of the genome
        chrom_length = chrom_features.length
        #If block method is selected
        if block:
            #null_sample is a list with each element as scores generated from a the block sampling method
            null_samples = block_null_sampling(chrom_features,points_process,chrom_length,win_size,int(chrom_length * block_coverage),int(1.0/block_coverage),num_sample)
        else:
            #Sampling from the NULL with random points as the point process
            null_samples = null_sampling(chrom_features,chrom_length,len(points_process),win_size,num_sample)
        #Sorts the null_sample such that it returns a list with each array as all the scores for a single point
        sorted_samples = sort_null_samples(null_samples)
        #computes the array forming the top and bottom percentage of the Null samples
        confidence_band = compute_confidence_band(sorted_samples,con_percent)
        #Takes the returned confidence band and seperated it to lower and upper
        upper_confidence_band = confidence_band[0]
        lower_confidence_band = confidence_band[1]
        #Applying differnt smoothing functions to the score and confidence band if indicated by user
        if aggregate_type:
            #If num_group was not specified
            if num_groups == 0:
                raise ValueError, "Need a nonzero value for number of groups if you wish to perform the moving_average aggregation"
            #Applies moving average to all three arrays-score, and two confidence bands
            if aggregate_type == 'moving_average':
                score_array = moving_average(score_array,num_groups)
                upper_confidence_band = moving_average(upper_confidence_band,num_groups)
                lower_confidence_band = moving_average(lower_confidence_band,num_groups)
            #Applies bin to all three arrays-score, and two confidence bands
            if aggregate_type == 'bin':
                score_array = bin(score_array,num_groups)
                upper_confidence_band = bin(upper_confidence_band,num_groups)
                lower_confidence_band = bin(lower_confidence_band,num_groups)
        
        #Plotting
        if plot_type == 'line':
            plot_score_confidence(score_array,upper_confidence_band,lower_confidence_band,num_sample)
        elif plot_type == 'bar':
            barplot_score_confidence(score_array,upper_confidence_band,lower_confidence_band,num_sample)
        elif plot_type == 'box':
            if aggregate_type == 'bin':
                #Allows for binning of the sorted_sample list such that the box plot will still work
                binned_list = []
                #nested for loops interate through each point in sorted sample
                for point_index in range(len(sorted_samples[0])):
                    sample = []
                    for sample_index in range(len(sorted_samples)):    
                        sample_point = sorted_samples[sample_index][point_index]
                        sample.append(sample_point)
                    #Bins the created sample
                    binned_list.append(bin(sample,num_groups))
                
                sorted_samples = sort_null_samples(binned_list)
                
            box_and_whisper_plot(score_array,sorted_samples)
        
        
        chrom_index += 1
    
    
def Inter2end(list_inter):
    """ Takes a feature in the region and returns a region object with
    only the end point of ther interval
    for ex: ([247, 357)] returns (357,357)"""
    
    #makes a copy  of the lis_iter
    new_list = list_inter[:]
    index = 0
    #Replaces the interval with the last value
    for intval in list_inter.iter_feature_regions():
        new_list[index] = base_types.interval(intval[1],intval[1])
        index = index + 1
    
    return new_list
        
def aggregate_scores(features,point_features,win_size,length):
    """The function computes the "scores" of the point process against the feature in genome A.
    this is done by returning ones for points along the window that are part of an interval in features
    and zeros where the point is not part of a feature.
    
    ex: if the point process was [5], the windowsize was 4 and features were [(1,3) (8,10)]
    then the new_array would be [1,1,0,0,0,0,1]. This process is repeated for each point and summed up
    to give a score
    
    The function's arguments are:
    -the feature
    -the point process
    -winsize and genome length
    
    The function returns the "score" generated """
    
    #Maxpoint is half the windowsize and the window stretches from -maxpoint to +maxpoing
    maxpoint = win_size/2
    #Initialize the final array to the zeros of the size of the window
    total_array = numpy.zeros(2*maxpoint+1)
    #Iterates through each point in the point process
    for point in point_features.iter_feature_regions():
        #sets the bounds on the window by adding and subtracting maxpoint from the point
        if point[0] != point[1]:
            raise ValueError, "For point_process file each interval must be the same value"
        rightpoint = point[0] + maxpoint
        leftpoint = point[0] - maxpoint
        #Forbid any points that when added with the window_size will be greater than the length of the genome
        if leftpoint < 0 or rightpoint > length:
            continue
        #subregion becomes a list with same length as window. any intervals present in features is retained
        #but shifted to the leftpoint = 0
        # gets the subregion- right index incremented by one to account for zero 
        subregion = features.get_subregion(leftpoint,rightpoint+1,True)            
        #To handle errors generated if there is no interval in subregion i.e adds (0,0) to a subregion of value []
        if not subregion:
            subregion = base_types.binary_region((),length = maxpoint*2+1)
            subregion.append((0,0))
            
        #Uses the bp_scored_array method in order to iterate through each point in the subregion putting ones where it is 
        #present in the feature and zero where it is not
        new_array = subregion.bp_score_array()

        #summing up the arrays
        total_array += new_array

    return total_array
        
def random_points_gen(genome_length,num_points):
    """Function accepts the genome length and num of points desired and 
    returns a binary_region object with randomly determined point process
    along the gnome"""
    import random
    #Creates an instances of a binary region with length = genome_length
    random_region = base_types.binary_region((),length = genome_length)
    #appends on a randon point along the genome for the num_points specified
    for num in range(num_points):
        # FIXME - changed to add
        random_point = int(random.randrange(0,genome_length))
        random_region.add(base_types.interval(random_point,random_point))
    return random_region


def null_sampling(features_A,chrom_length,num_points,win_size,num_sample):
    """
    This function takes sample "scores" from the null in order to compare the probability of
    attaining certain scores. This is accomplished by genereating random point arrays to use
    instead of the given point process B.
    
    The functions inputs are:
    -the features
    -genome_length
    -number of points, window size, and number of random runs used to test for confidence band
    
    The functon returns a list with each element as an array that contains the scores
    of sample generated from a random array."""
    #Initializing a list which will contain as each element the values generated from the compute_confidence_band
    #of random arrays
    scores_list = []
    for loop in range(num_sample):
        #This line appends on to the array_list a list of each scores from each run with a random points. 
        #We input the same values to aggregate_scores but for the point_features argument we send in 
        #points from the ran_points_gen function
        scores_list.append(aggregate_scores(features_A,random_points_gen(chrom_length,num_points),win_size,chrom_length))
    return scores_list
def block_null_sampling(features_A,points_B,chrom_length,win_size,block_size,num_block,nsamples):
    """This function performs sampling from the NULL using block techniques. It treats chunk of points from 
    the point process and computes the the score generated from moving the block to random place on 
    the features A tract
    
    this function accepts:
    -features from A
    -point process B
    -length of genome, window size
    -block size-how many basepairs does the block contain
    -number of samples to fun
    
    The function outputs:
    -a list with each elemeant an array containing scores generated from a sample
    """
    import random
    
    print 'using block_null_samples'
    sample_scores_list = []
    #creating nsamples number of sample
    for loop in range(nsamples):
        total_block_score = numpy.zeros(win_size/2 * 2 + 1)
        #For each block
        for loop in range(num_block):
            #Random place to create block            
            start_bp = random.randrange(0,chrom_length)
            #Gets Subregion
            block_region = points_B.get_subregion(start_bp,start_bp + block_size,True)
            #start point of the features where the block is going to move
            feature_start = random.randrange(0,chrom_length - block_size)
            #Shifting the points of the block to the start of the features that was randomly selected
            shift_points = base_types.binary_region((),length = chrom_length)
            for point in block_region.iter_feature_regions():
                shift_points.append((point[0] + feature_start, point[0] + feature_start))
            
            #Computes the score using aggregate scores
            total_block_score += aggregate_scores(features_A,shift_points,win_size,chrom_length)
       
        sample_scores_list.append(total_block_score)
    return sample_scores_list


def sort_null_samples(score_list):
    """Takes a list of generated scores from the null and sorts them such that each element of the list becomes 
    an array of all the scores for each individual point. and then proceeeds to sort them
    for ex: if the previous score was [[3,2,1],[4,6,4],[2,9,4]] the new list would be [[3,4,2],[2,6,9],[1,4,4]]"""
    final_list = []
    #Nested for loops to iterate through each score and each element of the score and assigns it to an element in the new array(final_list)
    for indexelement in range(len(score_list[0])):
        newlist = []
        for indexarray in range(len(score_list)):
            #For each list appends on to newlist the point in indexelement of that list
            newlist.append(score_list[indexarray][indexelement])
        #Attatches new_list which has the score for a particular point on the window to finali_list which 
        #will be a list of new_lists
        final_list.append(newlist)
        #Sorts each element inside final_list. uses list comprehension
        [points.sort() for points in final_list]

    return final_list
def compute_confidence_band(sorted_scores,percentage = .99):
    """This function takes a the sorted list generated form the funciton sort null samples and returns the confidence band as 
    determined by the percentage
    
    The function's inputs are:
    -a list of scores generated off the null that is then passed to sort_null_sampples
    -the percentage used to caluculate the confidence band
       ex: if percentage = .99 then the upper confident band will be greater than 99% of other values
       generated from it
    
    The function returns a list with two arrays
    -the first array contains the upper_confident_band
    -the second array contains the lower_confident_band
    """

    num_samples = len(sorted_scores[0])
    num_points =  len(sorted_scores)
    #Determines which point is the one such that it is greater than 99%of the other scores
    #Multiplies the number of runs by 1 percent. floors that value and adds 1 such that the index would be the 
    #correct one from the end
    #ex: if 100 runs then the index would be -(100 * .1 + 1) which is -2 and it's the second to last element.
    upper_confident_index = -int(numpy.floor((num_samples * (1-percentage)) + 1))
    lower_confident_index = -int(numpy.ceil(num_samples * percentage))
    #Initializes the bands
    upper_confident_band = numpy.zeros(num_points)
    lower_confident_band = numpy.zeros(num_points)
    count = 0
    #Sets the confident_band equal to the indexed value of the array in point_array
    for point_array in sorted_scores:
        upper_confident_band[count] = point_array[upper_confident_index]
        lower_confident_band[count] = point_array[lower_confident_index]
        count = count + 1
        
    return [upper_confident_band,lower_confident_band]


def moving_average(score_array,num_group):
    """
    Smoothing function that takes score_array and for each point averages it with the nearby values

    Function accepts an array and the number of groups to be used for the average. 
    example: if array was [3,2,1,4,3,2,4,6,3] and num_group was 4the first element of the returned array
    would be (3+2+1+4)/4
    
    Array is shortened in the process by a length of Num_group -1
    """

    #Sets the length of the new array generated. The num of group is /2 and multiplied by 2 such that 
    averaged_array =numpy.zeros(len(score_array)-(num_group -1))
    #Iterating through each of the averaged array
    for index in range(len(score_array )- (num_group -1)):
        #Setting each element of averaged array as the average of the score_array.
        averaged_array[index] = sum(score_array[index:index + num_group])/float(num_group)
    return averaged_array
        
def bin(score_array,num_group):
    """
    Smoothing function applied to an array of scores by summing groups of values in the array
    
    Accepts an array of scores and the number of groups
    
    ex: if array was [3,2,5,3,5,3] and the number of groups was 3 then the function will return an array [10,11]
    """
    #Size of the new array- is the original length divided by number of groups rounded down
    binned_array = numpy.zeros(len(score_array)/num_group)
    #Iterating through each index
    for index in range(len(score_array)/num_group):
        #The beginning index of where to start summing form the old array
        start_index = index * num_group
        #setting the element of binned_array equal to the sum of the new num_group of points in score_array
        binned_array[index] = sum(score_array[start_index:start_index + num_group])
    
    return binned_array
        
def barplot_score_confidence(score,upper,lower,num_sample):
    """
    Accepts a the scores generated formt he point process and the lower and upper confidence band
    
    Plots the score fromt the point process in a bargraph and creates confidence intervals 
    for each bar"""
    
    from rpy import r
    r.library('gplots')
    
    #Offering to save the plot or simply display it
    print "Enter a title for your plot if you wish to save it. Hit enter to just display the graph"
    title = sys.stdin.readline()[:-1]

    if title:
        #Saving file
        outfile = title
        r.bitmap(outfile,res = 200)
        r.barplot2(score,col = 'turquoise',ylab = 'Score', plot_ci = True, ci_u = upper, ci_l = lower)
        r.title(title)
        r.dev_off()
    else:
        r.barplot2(score,col = 'turquoise',ylab = 'Score', plot_ci = True, ci_u = upper, ci_l = lower)
        r.title('Barplot'+ '(Number of samples =' + str(num_sample) + ')')
        #Waiting for user response
        raw_input()
    

def plot_score_confidence(score,upper,lower,num_runs):
    """
    Accepts a the scores generated formt he point process and the lower and upper confidence band
    This function displays a plot of the score as a blue line and the confidence band as a red line"""
    from rpy import r
    
    combined = numpy.concatenate((score,upper,lower))
    #Determining the x axis of the plot. Shifted towards the left if the score is even number of values
    #due to smoothing it with aggregate types. Otherwise, incremented by one to account for the zero
    if len(score) % 2 == 1:
        x = range(-(len(score)/2),len(score)/2+1)
    else:
        x = range(-(len(score)/2),len(score)/2)
    #Offering to save the plot or simply display it
    print "Enter a title for your plot if you wish to save it. Hit enter to just display the graph"
    title = sys.stdin.readline()[:-1]

    if title:
        #Saving file
        outfile = title
        r.bitmap(outfile,res = 200)
        #Plots with limits equal to the high and low of the values present in the confidence and the score
        r.plot(x,score,col = 'blue',ylab = 'Score', type = 'b', ylim=(min(combined),max(combined)))
        r.lines(x, upper, col = 'red', lty='dashed',lwd = 3)
        r.lines(x, lower, col = 'red', lty='dashed',lwd = 3)
        r.legend("topright", legend=("Score","Confidence_Band", "Num of runs =" + str(num_runs)), col=('blue','red','green'),lty='dashed',lwd=3) 
        #Turning off device
        r.dev_off() 

    else:
        #Displaying the plot
        
        r.plot(x,score,col = 'blue',ylab = 'Score', type = 'b', ylim=(min(combined),max(combined)))
        #x = range(-25,26)
        r.lines(x, upper, col = 'red', lty='dashed')
        r.lines(x, lower, col = 'red', lty='dashed')
        r.legend("topright", legend=("Score","Confidence_Band"), col = ('blue','red'),lty='dashed',lwd=3) 
        r.title('Plot of score and confidence band' + "(Num of runs =" + str(num_runs)+ ')')
        #Waiting for user input to proceed
        raw_input()

def box_and_whisper_plot(score,sample_list):
    """generates a box and whisker plot the scores generated from the null at a given point"""
    from rpy import r

    # determine the minimum and maximum points for the purpose of setting the plot limits
    min_value = min(min(score), min( min(item) for item in sample_list ))
    max_value = max(max(score), max( max(item) for item in sample_list ))
    
    #Offering to save the plot or simply display it
    print "Enter a title for your plot if you wish to save it. Hit enter to just display the graph"
    title = sys.stdin.readline()[:-1]

    if title:
        #Saving file
        outfile = title
        r.bitmap(outfile,res = 200)
    
        # plot the NULL
        r.boxplot(sample_list,ylim = (min_value,max_value))
    
        # plot the real data
        r.lines(score,col = 'red', lty = 'dashed',lwd = 4)
        r.title(title)
        r.dev_off()
    else:
        # plot the NULL
        r.boxplot(sample_list,ylim = (min_value,max_value))
        # plot the real data
        r.lines(score,col = 'red', lty = 'dashed',lwd = 4)
        # wait for user input to continue
        r.title('boxplot')
        raw_input()
        

def usage():
    print >> output,    """
    THis program applies the block bootstrap method in order to determine 
    the significance of a point process relative to regions 
    on a second genome.
    
    
    The requied options to run this function are:
    -1, --feature_filename followed by the filename
    -2, --points_filename followed by the filename
    -l, --length_filename followed by the filename
    -w, --win_size to determine the window_size
    -r, --num_run for the the number of samples generated off the NUll
    
    Optional arguments are:
    -p, --plot_type to select type of plot: bar, line, box
    -c, --con_percent to determine the percentage applied to the confidence_band(ex:.99)
    -b, --block to use the block sampling method
    -s,, --block_size, Only needed if block is indicated: number to indicate the percentage of the genome covered by each block. 
    -a, --agg_type, To indicate the type of smoothing desired, moving_average or bin
    -g, --num_group, Number of points in each group for the smoothing

    
    for help type '-h', or '--help'(gives this message)
    
    To test with the Encode_annotations.bed 
    and segmentation.bed type -t or --test
    This will automatically perform the procedure with 
    conserved_sequences.bed ENCODE_annotations.bed, and 
    ENCODE_region_lengths.txt, a window of 50 10 samples, and 
    uses block sampling with .05 coverage Creates a boxplot
    
    """
def test(block= False):
    con_seq = open('conserved_sequences.bed','r')
    Encode_ann = open('ENCODE_annotations.bed','r')
    lengths_file = open('ENCODE_region_lengths.txt','r')

    features = base_types.parse_bed_file(con_seq,lengths_file)
    point_process =  base_types.parse_bed_file(Encode_ann,lengths_file)
    

    B_pointset = Inter2end(point_process.values()[0])
    compute_correlation_stat(features,B_pointset,50,10,.99,block=False,block_coverage = .05,plot_type = 'line',aggregate_type = 'moving_average',num_groups =4)

    

def main(argv):
    import getopt
    
    #Assigning default values to optional variables
    block = False
    block_size = 0
    aggregate_type = None
    num_groups = 0
    con_percent = .99
    plot_type = 'line'
    
    
    try:
        longargs = [ 'help','test','feature_filename=', 'points_filename', 'lengths_filename','win_size','block','num_runs=', 'plot=','con_percent','block_size','agg_type','num_groupss']
        opts,args = getopt.getopt(argv, 'ht1:2:l:w:br:p:c:s:a:g:',longargs)
     
    except getopt.GetoptError:
        usage()
        sys.exit()
    for opt,arg in opts:
        
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        if opt in ('-t','--test'):
            test()
            sys.exit()
        if opt in ('-1','--feature_filename='):
            features_filename = arg
        if opt in ('-2','--points_filename='):
            points_filename = arg
        if opt in ('-l','--lengths_filename='):
            lengths_filename = arg
        if opt in ('-w','--win_size'):
            win_size = int(arg)
        if opt in ('-r','--num_runs='):
            num_runs = int(arg)
        if opt in ('-p','--plot'):
            plot_type = arg
        if opt in ('-c', '--con_percent'):
            con_percent = float(arg)
        if opt in ('-b', '--block'):
            block = True
        if opt in ('-s', '--block_size'):
            block_size = float(arg)
        if opt in ('-a', '--agg_type'):
            aggregate_type = arg
        if opt in ('-g','--num_groups'):
            num_groups = int(arg)
    
    #Openining the files
    try:
        features_file = open(features_filename,'r')
        points_file = open(points_filename,'r')
        lengths_file = open(lengths_filename,'r')
    #If files had problems being opened
    except IOError:
        print 'The files provided could not be opened. Make sure you typed the correct name for each file'
        print 'type -h or --help for help on arguments for this program'
        sys.exit()
    #If not al the files were provided
    except UnboundLocalError:
        print 'Not all the files necesary were provided. Please type -h or --help for programs inputs'
        sys.exit()
    #Parsing files
    features = base_types.parse_bed_file(features_file,lengths_file)
    point_process =  base_types.parse_bed_file(points_file,lengths_file)
    
    #In case required information was not provided
    try:
        win_size, num_runs
    except UnboundLocalError:
        print 'Some of the required arguments were not inputed. type -h or --help for information on this programs arguments'
        sys.exit()
    #FIXME: deleted line if not using test data
    B_pointset = Inter2end(point_process.values()[0])
    #Sending to compute data
    compute_correlation_stat(features,B_pointset,win_size,num_runs,con_percent,block,block_size,plot_type,aggregate_type,num_groups)
        
                                
    
if __name__ == "__main__":
    main(sys.argv[1:])
