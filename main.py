import BlackLock_mod

def argparserLocal():
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='BlackLock', description='Visualize Genetic Sequence and Filter FASTQ and GZ file format')
    parser._optionals.title = 'optional arguments'
    
    subparsers = parser.add_subparsers(
        title='commands', description='Please choose command below:',
        dest='command'
    )
    subparsers.required = True

    visualize_command = subparsers.add_parser('VisualizeSeq', help='Generate HTML file for sequence summarizing')
    visualize_command.add_argument("-f", "--fname", type=str, default=None, help="Input Sequence filename (.fastq or .gz)")
    
    filter_command = subparsers.add_parser('FilterSeq', help='Filter FASTQ sequence and give output as a FASTQ file')
    filter_command.add_argument("-f", "--fname", type=str, default=None, help="Input Sequence filename (.fastq or .gz)") 
    filter_command.add_argument("-mi", "--minlen", type=int, default=0, help="Provide minimum length for filtering")
    filter_command.add_argument("-ma", "--maxlen", type=int, default=0, help="Provide maximum length for filtering")
    filter_command.add_argument("-ph", "--minPhred", type=int, default=0, help="Provide minimum Phred Score for filtering")
    filter_command.add_argument("-bar", "--barcode", type=str, default='a', help="Provide barcode name for filtering")
    filter_command.add_argument("-tb", "--totalb", type=int, default=0, help="Provide total base to be read for filtering")
    filter_command.add_argument("-seq", "--seq_search", type=str, default='', help="Provide sequence for searching for filtering")

    return parser

def main():
    parser = argparserLocal()
    args = parser.parse_args()
    #print(args)

    if args.command == 'VisualizeSeq':
        if args.fname == None:
            exit(parser.parse_args(['VisualizeSeq','-h']))
        BlackLock_mod.Omsin(args.fname)

    if args.command == 'FilterSeq':
        if args.fname == None:
            exit(parser.parse_args(['FilterSeq','-h']))
        BlackLock_mod.Tata2(args.fname, args.minlen, args.maxlen, args.minPhred, args.barcode, args.totalb, args.seq_search)