def make_samples_order(filein,order_name,fileout,tree=None):
    """
    Generates or appends to a samples-order.txt file for SAMPLES.db creation.

    This program appends orderings to a pre-existing samples-order.txt file. If
    the file does not exist it is created.

    PARAMETERS
    ----------
    filein : String
        filepath of the order of your samples (one sample per line or comma separated),
        or, the tree of your samples.
    tree : Boolean, default None
        state whether it a tree  
    """

    if tree == None:
        raise ValueError("Specify parameter `tree`.")

#   opening fileout for writing
    samples_order = open(fileout,"a+")

#   get some info from readonly
    readonly = open(fileout,"r").read()
    empty = True if len(readonly)==0 else False
    if order_name in readonly:
        raise ValueError("{} already exists as ordering".format(order_name))

#   read in order as string 
    ordering = open(filein,"r")
    ordering = ordering.read()

#   replaces newlines with commas (operation does not affect trees) or if 
#   format is already commas
    ordering = ordering.replace("\n",",")
    if ordering[-1]==",":
        ordering = ordering[:-1]

#   appends header if file is empty
    if empty:
        samples_order.write("attributes\tbasic\tnewick")

#   append samples_order info to fileout
    line = "\n"+order_name
    line += "\t\t"+ordering if tree else "\t"+ordering+"\t"
    samples_order.write(line)
