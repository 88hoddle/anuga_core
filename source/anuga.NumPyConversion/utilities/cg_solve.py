
import exceptions
class VectorShapeError(exceptions.Exception): pass
class ConvergenceError(exceptions.Exception): pass

##from numpy.oldnumeric import dot, array, Float, zeros
from numpy import dot, array, float, zeros
   
import logging, logging.config
logger = logging.getLogger('cg_solve')
logger.setLevel(logging.WARNING)

try:
    logging.config.fileConfig('log.ini')
except:
    pass

def conjugate_gradient(A,b,x0=None,imax=10000,tol=1.0e-8,iprint=0):
    """
    Try to solve linear equation Ax = b using
    conjugate gradient method

    If b is an array, solve it as if it was a set of vectors, solving each
    vector.
    """
    
    if x0 is None:
        x0 = zeros(b.shape, dtype=float)
    else:
        x0 = array(x0, dtype=float)

    b  = array(b, dtype=float)
    if len(b.shape) != 1 :
       
        for i in range(b.shape[1]):
            x0[:,i] = _conjugate_gradient(A, b[:,i], x0[:,i],
                                          imax, tol, iprint)
    else:
        x0 = _conjugate_gradient(A, b, x0, imax, tol, iprint)

    return x0
    
def _conjugate_gradient(A,b,x0=None,imax=10000,tol=1.0e-8,iprint=0):
   """
   Try to solve linear equation Ax = b using
   conjugate gradient method

   Input
   A: matrix or function which applies a matrix, assumed symmetric
      A can be either dense or sparse
   b: right hand side
   x0: inital guess (default the 0 vector)
   imax: max number of iterations
   tol: tolerance used for residual

   Output
   x: approximate solution
   """


   b  = array(b, dtype=float)
   if len(b.shape) != 1 :
      raise VectorShapeError, 'input vector should consist of only one column'

   if x0 is None:
      x0 = zeros(b.shape, dtype=float)
   else:
      x0 = array(x0, dtype=float)


   #FIXME: Should test using None
   if iprint == 0:
      iprint = imax

   i=1
   x = x0
   r = b - A*x
   d = r
   rTr = dot(r,r)
   rTr0 = rTr

   while (i<imax and rTr>tol**2*rTr0):
       q = A*d
       alpha = rTr/dot(d,q)
       x = x + alpha*d
       if i%50 :
           r = b - A*x
       else:
           r = r - alpha*q
       rTrOld = rTr
       rTr = dot(r,r)
       bt = rTr/rTrOld

       d = r + bt*d
       i = i+1
       if i%iprint == 0 :
          logger.info('i = %g rTr = %20.15e' %(i,rTr))

       if i==imax:
            logger.warning('max number of iterations attained')
            msg = 'Conjugate gradient solver did not converge: rTr==%20.15e' %rTr
            raise ConvergenceError, msg

   return x
