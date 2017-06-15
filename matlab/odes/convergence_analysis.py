#===============================================================================
# CLASS: Plain text database
#===============================================================================

class TextDB (object):
  """
  Plain text database, with data added line by line according to a user-defined 
  format. Name, description and format of each column is prescribed by using the
  SetField() method. When the new file is open, a header portion is added at the
  beginning of the file.
  
  Parameters
  ----------
  name : str
    Full file name.
  
  """
  def __init__(self,name):
    assert (type(name) is str)
    self._name     = name
    self._fields   = []
    self._ostream  = None
    self._template = None
  
  #-----------------------------------------------------------------------------
  def SetField (self, name, frmt, numtype=float, desc=None):
    """
    Add one field (i.e. one column) to the text database.
    
    Parameters
    ----------
    name : str
      Variable name (useful for debugging and extracting data later on).
    frmt : str
      Format (e.g. '+2.5e') for printing quantity to file.
    numtype : type
      Type of data.
    desc : str
      Description to be included in header section.
    
    """
    assert (type(name) is str)
    assert (type(frmt) is str)
    assert (type(desc) is str)  
    field = {'name': name, 'frmt': frmt, 'type': numtype, 'desc': desc }
    self._fields.append(field)
  
  #-----------------------------------------------------------------------------
  def open (self):
    """
    Open database stream. If opening for the first time, a new plain text file 
    is created with a header section at its beginning. If reopening, the 
    previous file is opened in 'append' mode.
    
    """
    # If opening for the first time, create a new file
    if self._ostream is None:
      self._ostream = open(self._name,'w')
      self._ostream.write (self._header())
      self._template = \
      ' '.join(['{:%s}'% field['frmt'] for field in self._fields]) + '\n'
    # If reopening, append to old file 
    else:
      self._ostream = open(self._ostream.name,'a')
  
  #-----------------------------------------------------------------------------
  def write (self,*data):
    """
    Add a new line to the database, according to the formatting defined in a 
    stored template.  The numerical values are passed as an argument list.
    
    """
    self._ostream.write(self._template.format(*data))
  
  #-----------------------------------------------------------------------------
  def close (self):
    """ Close database file stream.
    """
    self._ostream.close()
  
  #-----------------------------------------------------------------------------
  def _header (self):
    header = '# {:s}\n'.format(self._name)
    for idx,field in enumerate(self._fields):
      l0 = '#'
      l1 = '# Column {:2d}'
      l2 = '# ---------'
      l3 = '# {:s} {:s}'
      l4 = '# {:s}'
      lines = '\n'.join([l0,l1,l2,l3,l4]) + '\n'
      header += lines.format(idx,field['name'],field['type'],field['desc'])
    header += '\n'
    return header

#===============================================================================
# CLASS: Converter from Python float to LaTeX notation.
#===============================================================================

class Float2LatexConverter( object ):
  """
  Converter from Python number to string in LaTeX exponential notation.
  
  Number is first converted to a string in IEEE 754 floating-point format, 
  and then to a string in LaTeX-friendly exponential notation.
  
  Parameters
  ----------
  digits : int
    Number of digits to be kept in mantissa.
  
  Examples
  --------
  ::
    >>> convert = Float2LatexConverter( 5 )
    >>> convert( 1.25e-7 )
    '1.25000\\times 10^{-07}'
    >>> convert.set_prec( 1 )
    >>> convert( 1.25e-7 )
    '1.2\\times 10^{-07}'
  
  """
  def __init__( self, digits ):
    self.set_prec( digits )
  
  #-----------------------------------------------------------------------------
  def set_prec( self, digits ):
    """ Change number of digits in mantissa. """
    
    self._prec = digits
    self._str  = '{{:.{:d}e}}'.format( digits )
  
  #-----------------------------------------------------------------------------
  def __call__( self, x ):
    """ Convert number to LaTeX exponential notation. """
     
    float_str = self._str.format( x )
    return r'{0}\times 10^{{{1}}}'.format( *float_str.split('e') )

#===============================================================================
import math

#===============================================================================
# FUNCTION: Load data
#===============================================================================

def LoadData( *file_names ):
  """
  Load data from files.
  
  Returns
  -------
  mx : list of int
    Number of cells in domain for each of the simulations.
  
  data : list of dict
    For each file, dictionary of error norms for convergence analysis.
  
  """
  # Import numpy function for loading .dat files, and import ordered dictionary
  from numpy       import loadtxt
  from collections import OrderedDict
  
  # Load data and store it in a list of dictionaries
  data = []
  for fn in file_names:
    [mx, L1, L2, Li] = loadtxt( fn, unpack=True )
    record = OrderedDict()
    record['mx'] = mx.astype(int).tolist()
    record['L1'] = L1.tolist()
    record['L2'] = L2.tolist()
    record['Li'] = Li.tolist()
    data.append( record )

  # Check compatibility
  mx = data[0]['mx']
  for d in data:
    assert( d['mx'] == mx )
  
  return [mx, data]

#===============================================================================
# FUNCTION: Compute convergence ratios
#===============================================================================

def ComputeConvergence( mx, err ):
  """
  Compute convergence ratios.
  
  Parameters
  ----------
  mx : list of int
    Number of cells in domain for each of the simulations.
  
  err : list of floats
    Error norm for each of the simulations.
  
  """
  order = [None]
  for i in range(1,len(mx)):
    LogRatio_err = math.log(        err[i] /      err[i-1] )
    LogRatio_mx  = math.log( float(mx[i-1]) / float(mx[i]) )
    order.append( LogRatio_err / LogRatio_mx )
  
  return order

#===============================================================================
# FUNCTION: Create LaTeX table
#===============================================================================

def CreateLatexTable( mx, data,
                      file_name ='table.tex', 
                      mx_digits = 4,
                     err_digits = 2,
                     ord_digits = 2,
                     err_type   = 'L1'
                    ):
  
  # Number of runs and number of datasets
  NR = len(mx)
  ND = len(data)
  
  # Templates
  mx_fmt  = '{{:{0}d}}' .format(  mx_digits )
  err_fmt = '{:s}'
  ord_fmt = '{{:.{0}f}}'.format( ord_digits )
  
  row_template0 = ('${0}$' + ' & ${1}$ & ---' * ND + r'\\' + '\n') \
                  .format( mx_fmt, err_fmt )
  row_template  = ('${0}$' + ' & ${1}$ & ${2}$' * ND + r'\\' + '\n') \
                  .format( mx_fmt, err_fmt, ord_fmt )
  
  # Precomputed lines
  hline  = r'\hline'+'\n'
  begin  = r'\begin{tabular}{|r|'+'|c|c|'*ND+'}\n'
  end    = r'\end{tabular}'+'\n'
  header = r'\bf{Mesh} & ' + \
    ' & '.join( \
    [r'\bf{{M{} error}} & \bf{{Order}}'.format(i) for i in range(ND)]) + \
     r'\\' + '\n'
  
  # Exponential notation converter
  convert = Float2LatexConverter( err_digits )
  
  # Open output file
  f = open( file_name, 'w' )

  # Write table to file
  f.write( begin )
  f.write( hline )
  f.write( header)
  f.write( hline )
  f.write( hline )
  f.write( row_template0.format( mx[0], *[ convert(d[err_type][0]) for d in data ] ))
  f.write( hline )
  
  for i in range(1,NR):
    
    row_data = [ mx[i] ]
    for d in data:
      row_data += [ convert(d[err_type][i]), d['Order'][i] ]
    
    f.write( row_template.format( *row_data ) )
    f.write( hline )

  f.write( r'\end{tabular}' )
  f.write( '\n' )
  
  # Close output file
  f.close()

#===============================================================================
# FUNCTION: Parse input arguments
#===============================================================================

def parse_input():
  
  import argparse, sys
  
  parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description = '''Load error data from .dat files, estimate
                       convergence rates, and produce LaTeX table.''',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )
  
  parser.add_argument(metavar = 'FILE',
                      nargs   = '*',
                      dest    = 'files',
                      help    = 'input files')
  
  parser.add_argument('-o', '--output',
                      default = 'table.tex',
                      help    = 'output file name')
  
  parser.add_argument('-d', '--digits',
                      metavar = ('MX','ERR','ORDER'),
                      type    =  int,
                      nargs   =  3,
                      default = [4,2,2],
                      help    = "number of digits for 'mx', 'err', and 'order'")

  parser.add_argument('-e', '--err_type',
                      default = 'L1',
                      help    = """Error type for convergence Table.  Options
                      include L1, L2, and Li for the second through fourth
                      columns in the error tables.""")
  
  return parser.parse_args()

#===============================================================================
# FUNCTION: Main script
#===============================================================================

def main():
  
  # Parse input arguments
  args = parse_input()
  print(args)
  print('')
  
  # Load data
  mx, data = LoadData( *args.files )
  
  # Compute convergence ratios
  for d in data:
    d['Order'] = ComputeConvergence( mx, d[args.err_type] )
  
  # Create LaTeX table
  CreateLatexTable( mx, data,
                    file_name = args.output,
                    mx_digits = args.digits[0],
                   err_digits = args.digits[1],
                   ord_digits = args.digits[2])

#===============================================================================
if __name__ == '__main__':
  #Run as main program
  main()
