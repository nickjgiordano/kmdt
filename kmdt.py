# libraries
import re
from Bio import ExPASy, SwissProt, SeqIO
from tkinter import *

# function attached to button that submits an accession code
def accessionSearch():
    try:
        sInput = entryAccession.get() # get text field contents
        handle = ExPASy.get_sprot_raw(sInput) # for use in SwissProt.read method
        record = SwissProt.read(handle) # generates record from fasta code
        handle.close() # close handle since it's no longer in use
    except:
        # if exception is raised, display message to user
        lblResults.configure(text="invalid accession code!\n please try again...")
    else:
        # otherwise, submit sequence to motifFinder function
        motifFinder(record.sequence)
        
# function attached to button that submits a fasta code
def fastaSearch():
    sInput = entryFasta.get("1.0", END) # get text box contents
    # strip leading characters from fasta code
    if sInput.find("SV="):
        sInput = sInput[-(len(sInput)-sInput.find("SV=")-4):].strip()
    # remove line breaks from fasta code
    sequence = ""
    for line in sInput:
        sequence = sequence + line.strip()
    # submit code to motifFinder function
    motifFinder(sequence) 
    
# function that finds motif and displays it to user in a label
def motifFinder(sequence):
    # regex expression from brief, formatted for python
    searchPattern = "[^P][^PKRHW][VLSWFNQ][ILTYWFN][FIY][^PKRH]"
    # loops through sequence, creating matching motif substrings
    global match
    match = re.search(searchPattern, sequence)
    # variable to keep track of motif positions
    global lastFoundObject
    lastFoundObject = 1
    # default results string, before anything's added to it
    results = ""
    # if no matches are found, edit results string to display message
    if(match == None):
        results = "No results found!"
    # otherwise, loop through matches, adding motif sequences and positions to results string
    else:
        while(match != None):
            lastFoundObject += match.span()[0]
            results = results + match.group() + ", position {}\n".format(lastFoundObject)
            match = re.search(searchPattern, match.string[match.span()[1]:])
            lastFoundObject += 6;
    # add results string to label, removing trailing characters
    lblResults.configure( text=results.rstrip() )
    
# create gui window and set dimensions and position
root = Tk()
width = 600
height = 800
left = (root.winfo_screenwidth()-width)/2
top = (root.winfo_screenheight()-height)/2-20
root.geometry("%dx%d+%d+%d" % (width, height, left, top))
# colors and fonts
colorPurple = "#320064"
colorPurpleLight = "#5a2d87"
colorPurpleDark = "#280050"
fontPrimary = "Tahoma"
fontSecondary = "Rockwell"
fontMono = "Courier"
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
# set window title, background, and border
root.title("Keele Motif Detection Tool")
root.configure(bg=colorPurple, bd=12, relief="raised")
# create application title and info labels
Label(root, bg=colorPurple, font=(fontPrimary, 1) ).pack()
Label(root, bg=colorPurple, fg="white", font=(fontSecondary, 24), text="Keele Motif Detection Tool").pack()
Label(root, bg=colorPurple, font=(fontPrimary, 1) ).pack()
Label(root, bg=colorPurple, fg="white", font=(fontPrimary, 10), wraplength=540, text="Welcome! This is a tool designed especially for detecting protein regions that can form amyloid associations. The tool searches user-entered proteins for the following pattern motif:").pack()
Label(root, bg=colorPurple, font=(fontPrimary, 1) ).pack()
Label(root, bg=colorPurple, fg="white", font=(fontPrimary, 12), text="{P}1-{PKRHW}2-[VLSWFNQ]3-[ILTYWFN]4-[FIY]5-{PKRH}6".translate(subscript) ).pack()
Label(root, bg=colorPurple, font=(fontPrimary, 1) ).pack()
Label(root, bg=colorPurple, fg="white", font=(fontPrimary, 10), wraplength=540, text="Proteins can either be inputted via their UniProt accession code, or their FASTA code. Matching motifs will appear at the bottom of the window.").pack()
# create accession code text field
Label(root, bg=colorPurple, font=(fontSecondary, 20) ).pack()
Label(root, bg=colorPurple, fg="white", font=(fontSecondary, 14), text="Accession code").pack()
Label(root, bg=colorPurple, font=(fontSecondary, 1) ).pack()
entryAccession = Entry(root, bg=colorPurpleLight, fg="white", insertbackground="white", font=fontSecondary, width=20, justify="center", bd=4, relief="raised")
entryAccession.pack()
Label(root, bg=colorPurple, font=(fontSecondary, 1) ).pack()
btnAccession = Button(root, bg="black", fg="white", font=fontSecondary, width=20, relief="flat", text="Search", command=accessionSearch).pack()
# create fasta code text box
Label(root, bg=colorPurple, font=(fontSecondary, 20) ).pack()
Label(root, bg=colorPurple, fg="white", font=(fontSecondary, 14), text="FASTA code").pack()
Label(root, bg=colorPurple, font=(fontSecondary, 1)).pack()
entryFasta = Text(root, bg=colorPurpleLight, fg="white", insertbackground="white", font=(fontMono, 10), width=64, height=8, bd=4, relief="raised")
entryFasta.pack()
Label(root, bg=colorPurple, font=(fontSecondary, 1) ).pack()
btnFasta = Button(root, bg="black", fg="white", font=fontSecondary, width=57, relief="flat", text="Search", command=fastaSearch).pack()
# create results label
Label(root, bg=colorPurple, font=(fontSecondary, 20) ).pack()
Label(root, bg=colorPurple, fg="white", font=(fontSecondary, 14), text="Results").pack()
Label(root, bg=colorPurple, font=(fontSecondary, 1) ).pack()
lblResults = Label(root, bg=colorPurpleDark, fg="white", font=(fontMono, 9), width=40, bd=4, relief="raised", pady=10, text="enter a code to get results...")
lblResults.pack()
# ready window to be run
root.mainloop()
