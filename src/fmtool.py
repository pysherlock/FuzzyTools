#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
FMtool is a toolkit to manipulate fuzzy measures (aka capacities). It provides aggregation 
functions to perform a Multi-Criteria Decision Analysis (MCDA) on a set of criteria.

For more complex methods (such as Choquet and Sugeno), the FMtool class wraps a low-level implementation of an MCDA library 
written in plain C (see [1]), which is loaded in memory as a dynamic library (.dll, .dylib or .so) using ctypes. 

The toolkit also includes a function implementing the Ordered Weighted Average aggregation (OWA by Yager -- see [3]),
and two other methods that leverage the OWA operator to offer two generic aggregation functions: 
i) atLeast(x, k, n): aggregates a vector of criteria 'x' using a model that imposes that at least 'k' (out of 'n') criteria must be
satisfied (i.e., must be closer to 1)
ii) mostOf(x, n): aggregates a vector of criteria 'x' using a model that imposes that most of the 'n' criteria must be satisfied.

Example of use: see the test() function here under.
Written by: Olivier Thonnard 

Refs:
[1] G. Beliakov, A. Pradera, and T. Calvo. Aggregation Functions: A Guide for Practitioners.
Springer, Berlin, New York, 2007.
[2] V. Torra. The weighted OWA operator. Int. Journal of Intelligent Systems, 12(2):153–166, 1997.
[3] R. Yager. On ordered weighted averaging aggregation operators in multicriteria decision-making.
IEEE Trans. Syst. Man Cybern., 18(1):183–190, 1988.
[4] M. Grabisch and C. Labreuche. A decade of application of the choquet and sugeno integrals in multi-criteria decision aid. Annals of Operations Research, 2009.

"""

import numpy
import ctypes
import os, inspect
from scipy.interpolate import interp1d

BASE_PATH = os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0])


class FMtool:
    """This class provides methods to perform the aggregation of multiple criteria under a given decision-making scheme (MCDA).
	Some of methods are wrapping an implementation of a fuzzy measure toolkit written in plain C (given in Ref [1]).
	The fmtool class supports the calculation of OWA, Choquet and Sugeno integrals (Weighted OWA to come soon).
	Some methods are also provided to construct lambda-fuzzy and 2-additive fuzzy measures.
	"""
    alpha_RIM = {'atleast1': 0.001, 'few': 0.1, 'some': 0.5, 'half': 1, 'many': 2, 'most': 10, 'all': 1000}

    def __init__(self, nr_of_criteria, aggr_type='ChoquetInt'):
        self._n = nr_of_criteria
        self._m = 2 ** nr_of_criteria
        self._mp = numpy.asarray([self._m], dtype=numpy.uint)
        self._v = numpy.zeros(self._m, dtype=numpy.double)
        self._mobius = numpy.zeros(self._m, dtype=numpy.double)
        self._lambda_value = 0
        self._w = []
        self._p = []
        self._interactions = {}
        self._importance_factors = {}
        self._aggr_type = aggr_type

        # set default values for importance factors (ie, Shapley values)
        for k in range(1, nr_of_criteria + 1):
            self._importance_factors[k] = 0.0

        # set default values for interaction indices (ie, no interaction = 0)
        for i in range(1, nr_of_criteria):
            for j in range(i + 1, nr_of_criteria + 1):
                self._interactions[(i, j)] = 0.0

                # #try:
                # #	self._libfmtools = numpy.ctypeslib.load_library('libfuzzymeasuretools', BASE_PATH)
                # #except OSError:
                # #	self._libfmtools = numpy.ctypeslib.load_library('libfuzzymeasuretools', os.getcwd())
                #
                # # Initialization of the shared library is needed!
                # # C++ implementation: void Preparations_FM(int n, unsigned int *m)
                # self._libfmtools.Preparations_FM.argtypes = [ctypes.c_int, numpy.ctypeslib.ndpointer(dtype = numpy.uint)]
                # self._libfmtools.Preparations_FM.restype = ctypes.c_void_p
                # self._libfmtools.Preparations_FM(self._n,self._mp)
                #
                # # We must set the correct data types for the functions arguments and the return variables
                # # Choquet C++ implementation: double Choquet(double*x, double* v, int n, unsigned int m)
                # self._libfmtools.Choquet.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double),numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.Choquet.restype = ctypes.c_double
                #
                # # Mobius C++ implementation: void Mobius(double* v, double* Mob, int n,unsigned int m)
                # self._libfmtools.Mobius.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double),numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.Mobius.restype = ctypes.c_void_p
                #
                # # ChoquetMob C++ implementation: double ChoquetMob(double*x, double* Mob, int n, unsigned int m)
                # self._libfmtools.ChoquetMob.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double),numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.ChoquetMob.restype = ctypes.c_double
                #
                # # void ConstructLambdaMeasure(double *singletons, double *lambda, double *v, int n,unsigned int m)
                # self._libfmtools.ConstructLambdaMeasure.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), numpy.ctypeslib.ndpointer(dtype = numpy.double), numpy.ctypeslib.ndpointer(dtype = numpy.double),ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.ConstructLambdaMeasure.restype = ctypes.c_void_p
                #
                # # Sugeno C++ implementation: double Sugeno(double*x, double* v, int n, unsigned int m)
                # self._libfmtools.Sugeno.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double),numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.Sugeno.restype = ctypes.c_double
                #
                # # Orness C++ implementation: double Orness(double* Mob,  int n, unsigned int m)
                # self._libfmtools.Orness.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.Orness.restype = ctypes.c_double
                #
                # # Entropy C++ implementation: double Entropy(double* v,  int n, unsigned int m)
                # self._libfmtools.Entropy.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.Entropy.restype = ctypes.c_double
                #
                # # Shapley C++ implementation: void Shapley(double* v, double* x, int n, unsigned int m)
                # self._libfmtools.Shapley.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double),numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.Shapley.restype = ctypes.c_void_p
                #
                # # Banzhaf C++ implementation: void Banzhaf(double* v, double* x, int n, unsigned int m)
                # self._libfmtools.Banzhaf.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double),numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.Banzhaf.restype = ctypes.c_void_p
                #
                # # Interaction C++ implementation: void Interaction(double* Mob, double* w, unsigned int m)
                # self._libfmtools.Interaction.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double),numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_uint]
                # self._libfmtools.Interaction.restype = ctypes.c_void_p
                #
                # # Defining all functions: isMeasure*
                # # C++ implementation: int IsMeasureBalanced(double* v,unsigned int m)
                # self._libfmtools.IsMeasureBalanced.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_uint]
                # self._libfmtools.IsMeasureBalanced.restype = ctypes.c_int
                #
                # # C++ implementation: int IsMeasuresubadditive(double* v,unsigned int m)
                # self._libfmtools.IsMeasuresubadditive.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_uint]
                # self._libfmtools.IsMeasuresubadditive.restype = ctypes.c_int
                #
                # # C++ implementation: IsMeasuresuperadditive(double* v,unsigned int m)
                # self._libfmtools.IsMeasuresuperadditive.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_uint]
                # self._libfmtools.IsMeasuresuperadditive.restype = ctypes.c_int
                #
                # # C++ implementation: int IsMeasuresupermodular(double* v,unsigned int m)
                # self._libfmtools.IsMeasuresupermodular.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_uint]
                # self._libfmtools.IsMeasuresupermodular.restype = ctypes.c_int
                #
                # # C++ implementation: int IsMeasuresubmodular(double* v,unsigned int m);
                # self._libfmtools.IsMeasuresubmodular.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_uint]
                # self._libfmtools.IsMeasuresubmodular.restype = ctypes.c_int
                #
                # # C++ implementation: int IsMeasureadditive(double* v, int n,unsigned int m)
                # self._libfmtools.IsMeasureadditive.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.IsMeasureadditive.restype = ctypes.c_int
                #
                # # C++ implementation: int IsMeasuresymmetric(double* v, int n,unsigned int m)
                # self._libfmtools.IsMeasuresymmetric.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_int, ctypes.c_uint]
                # self._libfmtools.IsMeasuresymmetric.restype = ctypes.c_int
                #
                # # C++ implementation: int IsMeasureselfdual(double* v,unsigned int m)
                # self._libfmtools.IsMeasureselfdual.argtypes = [numpy.ctypeslib.ndpointer(dtype = numpy.double), ctypes.c_uint]
                # self._libfmtools.IsMeasureselfdual.restype = ctypes.c_int

    #
    def getAggrType(self):
        return self._aggr_type

    def setAggrType(self, aggr_type):
        self._aggr_type = aggr_type

    def setParams(self, params):
        """
		Convenience method to set parameters of the FMTool object for the aggregation method defined during initialization.
		The params argument must be a dictionary containing all required parameters associated to the given aggregation method.
		"""
        if self._aggr_type == 'OWA':
            w = params.get('w')
            if type(w) in [list, numpy.ndarray]:
                self.setOWAweights(w)
            if type(w) in [str, unicode]:
                if w.lower() in self.alpha_RIM.keys():
                    self.setRIMQuantifier(w)
                elif w.lower() == 'atleastk':
                    k = int(params.get('k'))
                    self.setAtLeastWeights(k)
        elif self._aggr_type == 'WOWA':
            w = params.get('w')
            if type(w) in [list, numpy.ndarray]:
                self.setOWAweights(w)
            if type(w) in [str, unicode]:
                if w.lower() in self.alpha_RIM.keys():
                    # print "setting alpha RIM: %s"%w.lower() # debug
                    self.setRIMQuantifier(w)
                elif w.lower() == 'atleastk':
                    k = int(params.get('k'))
                    # print "setting 'atleastk' vector with k=%d"%k # debug
                    self.setAtLeastWeights(k)
            p = params.get('p')
            # print "setting WOWA vector: %s"%str(p) # debug
            self.setWOWAWeights(p)
        elif self._aggr_type == 'Choquet':
            v = params.get('v')
            self.setFuzzyMeasure(v)
        elif self._aggr_type == 'ChoquetInt':
            importance_factors = params.get('importance_factors')
            synergies = params.get('synergies')
            redundancies = params.get('redundancies')
            for key, val in importance_factors.items():
                self.setImportance(int(key), val)
            if redundancies and type(redundancies) == dict:
                for key, val in redundancies.items():
                    key = eval(key)  # to re-transform the string '(x,y)' into a Python tuple
                    self.setInteraction(key, abs(val), 'redund')
            if synergies and type(synergies) == dict:
                for key, val in synergies.items():
                    key = eval(key)
                    self.setInteraction(key, val, 'syn')
        elif self._aggr_type == 'ChoquetInt2':
            importance_factors = params.get('importance_factors')
            interactions = params.get('interactions')
            for key, val in importance_factors.items():
                self.setImportance(int(key), val)
            for key, val in interactions.items():
                key = eval(key)
                self.setRawInteraction(key, val)
        else:
            self.setMostOfWeights()

    def setImportance(self, crit, value):
        """
		Sets an importance factor to a criteria. The value must be in [0,1].
		"""
        if value < 0 or value > 1:
            print "error: an importance factor must be in the range [0,1]"
            return
        if crit not in range(1, self._n + 1):
            print "error: criterion must be in range [1,%d]" % self._n

        self._importance_factors[crit] = value

        # load balance the other factors if the new importance factor is too high
        if sum(self._importance_factors.values()) > 1.0:
            excess = (sum(self._importance_factors.values()) - 1.0) / (self._n - 1)
            for key, val in self._importance_factors.items():
                if key != crit:
                    self._importance_factors[key] = val - excess

    def setRawInteraction(self, pair, value):
        """
		Sets a raw interaction value for a pair of criteria (I_ij). Does not check for consistency of these values.
		"""
        self._interactions[pair] = value

    def setInteraction(self, pair, ratio, interaction_type='syn'):
        """
		Assigns an interaction ratio to a pair of criteria. The 'interaction_type' can be either 'syn' ( for synergy, or positive interaction),
		or 'redund' (for redundancy, or negative interaction). The 'ratio' (a real number in [0,1]) indicates how much of interaction we want to assign.
		Note also that this method will check for contraints with respect to all other interactions involving the same criteria, 
		and, if needed, it will load balance with other interaction indices to comply with the global monotonicity requirement for the capacity.
		"""
        if len(pair) != 2 or not isinstance(pair, tuple):
            print 'error: a pair of criteria must be of the form: (i,j) [for example: (1,2)]'

        crit1 = pair[0]
        crit2 = pair[1]
        # verify first what is the max. range that we can set for all interaction indices of these criteria
        max_range = min(2 * self._importance_factors[crit1], 2 * self._importance_factors[crit2])
        interaction = max_range * abs(ratio)
        if interaction_type == 'redund':
            self._interactions[pair] = -interaction
        else:
            self._interactions[pair] = interaction

        # we must verify that we don't violate the constraints regarding the other interactions involving those criteria
        sum_crit1 = 0.0
        other_crit1 = []
        sum_crit2 = 0.0
        other_crit2 = []

        for key, val in self._interactions.items():
            if crit1 in key and crit2 not in key:
                sum_crit1 += abs(val)
                other_crit1.append(key)
            if crit2 in key and crit1 not in key:
                sum_crit2 += abs(val)
                other_crit2.append(key)
        # print 'crit %d: %g, %s' % (crit1, sum_crit1, str(other_crit1))
        # print 'crit %d: %g, %s' % (crit2, sum_crit2, str(other_crit2))
        # check constraints for criterion 1
        if (interaction > 2 * self._importance_factors[crit1] - sum_crit1):
            excess = abs(2 * self._importance_factors[crit1] - sum_crit1 - interaction)
            # print 'warning: must compensate I%d%d for criterion %d, excess: %g' % (crit1,crit2, crit1, excess)
            for k in other_crit1:
                if self._interactions[k] == 0.0: continue
                if abs(self._interactions[k]) >= excess:
                    if self._interactions[k] > 0:
                        self._interactions[k] -= excess
                    else:
                        self._interactions[k] += excess
                    break
                else:
                    self._interactions[k] = 0.0
                    excess -= abs(self._interactions[k])
        # check constraints for criterion 2
        if (interaction > 2 * self._importance_factors[crit2] - sum_crit2):
            excess = abs(2 * self._importance_factors[crit2] - sum_crit2 - interaction)
            print 'warning: must compensate I%d%d for criterion %d, excess: %g' % (crit1, crit2, crit2, excess)
            for k in other_crit2:
                if self._interactions[k] == 0.0: continue
                if abs(self._interactions[k]) >= excess:
                    if self._interactions[k] > 0:
                        self._interactions[k] -= excess
                    else:
                        self._interactions[k] += excess
                    break
                else:
                    self._interactions[k] = 0.0
                    excess -= abs(self._interactions[k])
                    # return self._interactions

    def Aggregate(self, z):
        """
		Returns the aggregated value of input vector z. 
		The aggregation method was defined during the initialization of the FMtool object or using the method 'setAggrType'.
		"""
        if self._aggr_type == 'OWA':
            return self.OWA(z)
        elif self._aggr_type == 'WOWA':
            return self.WOWA(z)
        elif self._aggr_type == 'ChoquetInt':
            return self.ChoquetInt(z)
        elif self._aggr_type == 'Choquet':
            return self.Choquet(z)
        elif self._aggr_type == 'ChoquetInt2':
            return self.ChoquetInt(z)
        else:
            return self.OWA(z)

    def ChoquetInt(self, z, verbose=False):
        """
		Computes the Choquet integral based on following formula: 
		C_v(z) = \sum_{i,j \in N \vert I_{ij} > 0} ( z_i \wedge z_j) I_{ij} 
		       + \sum_{i,j \in N \vert I_{ij} < 0} ( z_i \vee z_j) \vert I_{ij} \vert 
		       + \sum_{i \in N} z_i \left[ \phi_{i} - \frac{1}{2} \sum_{j \neq i } \vert I_{ij} \vert \right]
		This is formula (2.6) given by grabisch in Ref [4].
		"""
        synergies = {}
        redundancies = {}
        importances = {}
        sum_Iij = {}

        for i in range(1, self._n + 1):
            sum_Iij[i] = 0.0

        for key, val in self._interactions.items():
            if val == 0.0: continue

            if val > 0.0:
                contrib = min(z[key[0] - 1], z[key[1] - 1]) * val
                if contrib > 0.0: synergies[key] = contrib
                for i in range(1, self._n + 1):
                    if i in key: sum_Iij[i] += val

            if val < 0.0:
                contrib = max(z[key[0] - 1], z[key[1] - 1]) * abs(val)
                if contrib > 0.0: redundancies[key] = contrib
                for i in range(1, self._n + 1):
                    if i in key: sum_Iij[i] += abs(val)

        for i in range(1, self._n + 1):
            if z[i - 1] == 0.0:
                continue
            else:
                contrib = z[i - 1] * (self._importance_factors[i] - 0.5 * sum_Iij[i])
                if contrib > 0.0: importances[i] = contrib

        choq = sum(synergies.values()) + sum(redundancies.values()) + sum(importances.values())

        if verbose:
            return choq, synergies, redundancies, importances
        else:
            return choq

    def setFuzzyMeasure(self, fuzzy_measure):
        """
		Sets the fuzzy measure v. By convention, each element of the 2^n input vector represents the importance of a subset (or coalition) of criteria, in the following order:
		--> (empty set), (1), (2), (1,2), (3), (1,3), (2,3), (1,2,3), (4), (1,4), (2,4), (1,2,4), (3,4), (1,3,4), (2,3,4), (1,2,3,4), etc.
		(this convention comes from the internal binary representation of all those combinations, ie for n=4:
		(0000), (0001), (0010), (0011), (0100), (0101), (0110), (0111), ... , (1111), etc)
		"""
        self._v = numpy.asarray(fuzzy_measure, dtype=numpy.double)
        # already compute the Mobius representation of 'v' as soon it is defined.
        self._Mobius()

    def getFuzzyMeasure(self):
        return self._v

    def printFuzzyMeasure(self):
        """ Convenience function that prints a formatted view of the defined fuzzy measure. """
        pass

    def getLambda(self):
        return self._lambda_value

    def constructLambdaMeasure(self, singletons):
        """
		Given the values of the fuzzy measure at singletons, finds the appropriate
		lambda, and constructs the rest of the fuzzy measure. Returns lambda and v as output
		"""
        # void ConstructLambdaMeasure(double *singletons, double *lambda, double *v, int n,unsigned int m)
        singletons = numpy.asarray(singletons, dtype=numpy.double)
        lambda_p = numpy.asarray([self._lambda_value], dtype=numpy.double)
        self._libfmtools.ConstructLambdaMeasure(singletons, lambda_p, self._v, self._n, self._m)
        self._lambda_value = lambda_p[0]
        self._Mobius()

        return (self._v, self._lambda_value)

    def Choquet(self, x):
        """Calculates the (discrete) Choquet integral of the vector 'x', wrt the fuzzy measure v (defined as class member)."""

        data = numpy.asarray(x, dtype=numpy.double)
        choq = self._libfmtools.Choquet(data, self._v, self._n, self._m)

        return choq

    def Sugeno(self, x):
        """Calculates the value of a discrete Sugeno integral of vector 'x', wrt fuzzy measure v (defined as class member)."""

        data = numpy.asarray(x, dtype=numpy.double)
        sug = self._libfmtools.Sugeno(data, self._v, self._n, self._m)

        return sug

    def ChoquetMob(self, x):
        """ This is an alternative calculation of the Choquet integral using the Mobius transform.
		    It is not as efficient as Choquet(self, x). Provided for testing purposes. """

        data = numpy.asarray(x, dtype=numpy.double)
        choq = self._libfmtools.ChoquetMob(data, self.getMobius(), self._n, self._m)

        return choq

    def Orness(self):
        """Calculates orness value of the fuzzy measure v (defined as class member)."""

        # double Orness(double* Mob,  int n,unsigned int m)
        return self._libfmtools.Orness(self.getMobius(), self._n, self._m)

    def Entropy(self):
        """Calculates the entropy of the fuzzy measure v (defined as class member)."""

        # double Entropy(double* v,  int n,unsigned int m)
        return self._libfmtools.Entropy(self._v, self._n, self._m)

    def getMobius(self):
        """Returns the Mobius representation of the fuzzy measure."""

        return self._mobius

    def _Mobius(self):
        """
		Calculates the Mobius representation of the fuzzy measure v (defined as class member).
		"""
        self._libfmtools.Mobius(self._v, self._mobius, self._n, self._m)
        return None

    def Shapley(self):
        """Calculates and returns the array of Shapley values for the defined fuzzy measure v."""

        # void Shapley(double* v, double* x, int n,unsigned int m)
        shap = numpy.zeros(self._n, dtype=numpy.double)
        self._libfmtools.Shapley(self._v, shap, self._n, self._m)

        return shap

    def Banzhaf(self):
        """Calculates and returns the array of Banzhaf values for the defined fuzzy measure v."""

        # void Banzhaf(double* v, double* x, int n,unsigned int m)
        banz = numpy.zeros(self._n, dtype=numpy.double)
        self._libfmtools.Banzhaf(self._v, banz, self._n, self._m)

        return banz

    def Interaction_from_Mobius(self):
        """Returns all 2^n interaction indices of the defined fuzzy measure v."""

        # void Interaction(double* Mob, double* w,unsigned int m)
        interactions = numpy.zeros(self._m, dtype=numpy.double)
        self._libfmtools.Interaction(self.getMobius(), interactions, self._m)

        return interactions

    def isMeasureBalanced(self):
        # int IsMeasureBalanced(double* v,unsigned int m)
        return self._libfmtools.IsMeasureBalanced(self._v, self._m)

    def isMeasuresubadditive(self):
        # int IsMeasuresubadditive(double* v,unsigned int m)
        return self._libfmtools.IsMeasuresubadditive(self._v, self._m)

    def isMeasuresuperadditive(self):
        # int IsMeasuresuperadditive(double* v,unsigned int m)
        return self._libfmtools.IsMeasuresuperadditive(self._v, self._m)

    def isMeasuresupermodular(self):
        # int IsMeasuresupermodular(double* v,unsigned int m)
        return self._libfmtools.IsMeasuresupermodular(self._v, self._m)

    def isMeasuresubmodular(self):
        # int IsMeasuresubmodular(double* v,unsigned int m)
        return self._libfmtools.IsMeasuresubmodular(self._v, self._m)

    def isMeasureadditive(self):
        # int IsMeasureadditive(double* v, int n,unsigned int m)
        return self._libfmtools.IsMeasureadditive(self._v, self._n, self._m)

    def isMeasuresymmetric(self):
        # int IsMeasuresymmetric(double* v, int n,unsigned int m)
        return self._libfmtools.IsMeasuresymmetric(self._v, self._n, self._m)

    def isMeasureselfdual(self):
        # int IsMeasureselfdual(double* v,unsigned int m)
        return self._libfmtools.IsMeasureselfdual(self._v, self._m)

    def setAtLeastWeights(self, k):
        """
		Calculates and sets the appropriate weighting vector for the OWA aggregation by applying the rule 'at least k (out of n) criteria'.
		In practise, this means that the top 'k' criteria (of the ordered input x) will be given weights according to a normal distribution centered on mu=k, 
		with sigma=1 and sum(w_i) ~ 0.69 for i in [1:k].
		This is a convenience function that allows one to let an fmtool object set the weighting only once, and then reuse it for many computations.
		So, this is exactly the same as what happens in the atLEast method, except that here we don't return the OWA aggregated value.
		"""

        mu = k
        sigma = 1

        X = numpy.arange(1, self._n + 1)
        W = numpy.exp(-0.5 * ((X - mu) / sigma) ** 2) / (numpy.sqrt(2 * numpy.pi) * sigma)
        residual = 1.0 - sum(W)
        if residual > 0:
            W[k - 1] = W[k - 1] + residual
        self._w = W

    def setOWAweights(self, w):
        """
		This function simply sets the weights used in (W)OWA aggregation methods.
		"""
        self._w = numpy.asarray(w, dtype=numpy.float32)

    def setRIMQuantifier(self, kind='half'):
        """
		Adjusts the OWA weights to a given fuzzy linguistic schema: 
		- atleast1
		- few
		- some
		- half (default)
		- many
		- most
		- all
		"""
        schema = kind.lower()
        alpha = self.alpha_RIM.get(schema)
        if alpha:
            N = self._n
            w = numpy.zeros(N)

            for j in range(1, N + 1):
                w[j - 1] = (float(j) / N) ** alpha - ((float(j - 1)) / N) ** alpha

            self.setOWAweights(w)

    def setWMWeights(self, p):
        """
		Sets the weighting vector 'p' used in WM.
		"""
        self._p = numpy.asarray(p, dtype=numpy.float32)

    def helpDefineImportances(self, **kwargs):
        """
		Helper method to automatically define the importance weights as used in WOWA (the 'p' vector).
		Syntax for keyword arguments is: c1='h', c2='hh', c3='l', ... (where 'hh'=very high, 'h'=high, 'l'=low, 'll'=very low)
		"""
        avg_val = 1.0 / self._n
        p = [avg_val for i in range(0, self._n)]
        crit_list = ['c' + str(i) for i in range(1, self._n + 1)]

        if len(kwargs) > 0:
            for k, v in kwargs.items():
                if k in crit_list:

                    idx = crit_list.index(k)

                    if v.lower() == 'h':
                        p[idx] = avg_val * 1.5

                    if v.lower() == 'hh':
                        p[idx] = avg_val * 2.0

                    if v.lower() == 'l':
                        p[idx] = avg_val / 2.0

                    if v.lower() == 'll':
                        p[idx] = avg_val / 3.0

        tot_weight = sum(p)
        for i, v in enumerate(p):
            p[i] = v / tot_weight

        return p

    def setWOWAWeights(self, p):
        """
		Sets the weighting vector 'p' used in WOWA. The other weighting vector 'w' needs to be set first, as these two are used
		together to calculate w* (the interpolation function), which is needed for the determination of the final weights.
		"""
        x_data = [0]
        y_data = [0]
        self._p = numpy.asarray(p, dtype=numpy.float32)
        # self._psigma = numpy.array(sorted(p, reverse=True))
        w = self._w
        N = self._n
        for i in range(1, N + 1):
            x_data.append(float(i) / N)
            y_data.append(sum(w[:i]))

        self._xdata = x_data
        self._ydata = y_data
        w_star = interp1d(x_data, y_data)
        self.w_star = w_star

    def setMostOfWeights(self):
        """
		Calculates and sets the appropriate weighting vector for the OWA aggregation by applying the rule 'at least k (out of n) criteria'.
		"""

        k = numpy.round(0.8 * self._n)
        if k < 1.0:
            k = 1.0
        self.setAtLeastWeights(k)

    def mostOf(self, x):
        """
		Calculates the aggregation of the criteria given an input vector x by applying the rule 'most of the n criteria must be satisfied'.
		In practise, this means that we take 80%  of the criteria into account to find the last relevant input of the input vector x (denoted by k). 
		For example, with n=5 we will find that k = 0.80 * 5 = 4. Then, we apply the function atLeast(x, k, n) with the value k just found (see the doc of that function for more details).
		"""

        k = numpy.round(0.8 * self._n)
        if k < 1.0:
            k = 1.0
        return self.atLeast(x, k)

    def atLeast(self, x, k):
        """
		Calculates the aggregation of the criteria given in input vector x by applying the rule 'at least k (out of n) criteria must be satisfied'.
		In practise, this means that the top 'k' criteria (of the ordered input x) will be given weights according to a normal distribution centered on mu=k, 
		with sigma=1 and sum(w_i) ~ 0.69 for i in [1:k].	
		"""
        mu = k
        sigma = 1

        X = numpy.arange(1, self._n + 1)
        W = numpy.exp(-0.5 * ((X - mu) / sigma) ** 2) / (numpy.sqrt(2 * numpy.pi) * sigma)
        self._w = W

        return self.OWA(x)

    def WM(self, x):
        """
		Calculates the Weighted Mean of vector 'x' using a predefined weighting vector 'p'.
		"""
        x = numpy.array(x)

        return numpy.dot(x, self._p)

    def OWA(self, x):
        """
		Calculates the Ordered Weigted Average of vector 'x' using a predefined weighting vector 'w', as presented by Yager (See Ref.[3])
		"""
        xs = numpy.array(sorted(x, reverse=True))

        return numpy.dot(xs, self._w), self._w

    def WOWA(self, x):
        """
		Calculates the Weighted OWA aggregation of 'x' using weighting vectors 'w' and 'p' (according to Torra's implementation, see Ref.[2])
		"""
        xs = numpy.array(sorted(x, reverse=True))
        values_cache = []
        sigma = []

        for j in xs:
            if j in values_cache:
                continue
            idx = numpy.nonzero(x == j)[0]
            sigma.extend(idx)
            values_cache.append(j)

        print sigma, values_cache
        # print "sigma: %s"%str(sigma)
        N = self._n
        p = self._p
        w_star = self.w_star
        # print N, p, w_star
        omega = numpy.zeros(N)

        omega[0] = w_star(p[sigma[0]])

        for i in range(2, N + 1):
            part1 = min(sum(p[sigma[:i]]), 1)
            part2 = min(sum(p[sigma[:i - 1]]), 1)
            omega[i - 1] = w_star(part1) - w_star(part2)

        return numpy.dot(xs, omega), omega


def binomialCoefficient(n, k):
    if k < 0 or k > n:
        return 0
    if k > n - k:  # take advantage of symmetry
        k = n - k
    c = 1
    for i in range(k):
        c = c * (n - (k - (i + 1)))
        c = c // (i + 1)
    return c

def test_WOWA():
    weight = [2, 1.5, 1, 0.8, 0.7, 0.5, 0.3, 0.2]
    Sum = sum(weight)
    weight = [weight[i] / Sum for i in range(len(weight))]

    fmt = FMtool(8, 'Choquet')
    fmt.setRIMQuantifier('atleast1')
    fmt.setOWAweights([0.25702141884872337, 0.42375668028141705, 0.25702141884872337, 0.05734923034929231,
                       0.004707511494303314, 1.4215452955844868E-4, 1.5791941767602698E-6, 6.453805657372647E-9])
    # fmt.setWOWAWeights([2, 1.5, 1, 0.8, 0.7, 0.5, 0.3, 0.2])
    fmt.setWOWAWeights(weight)

    print "WOWA test 1: ", fmt.WOWA([0.982, 0, 0.6224593312018546, 0.5740969929676947, 1, 0.799, 0, 0.819])
    print "WOWA test 2: ", fmt.WOWA([0.984, 0, 0, 0.014926521817160062, 0, 0.44, 0, 0.657])
    print "WOWA test 3: ", fmt.WOWA([0.824, 0, 0, 0, 0.07, 0.108, 0, 0.696])
    print "WOWA test 4: ", fmt.WOWA([0.753, 0, 0.6224593312018546, 0.5740969929676947, 1, 0.799, 0, 0.819])
    print "WOWA test 5: ", fmt.WOWA([0.853, 0, 0, 0, 1, 0.454, 0, 1]) # 0.975

    print "WOWA test 6: ", fmt.WOWA([0.408, 0, 0.6224593312018546, 0, 0, 0.475, 0, 1]) #0.534:0.578
    print "WOWA test 7: ", fmt.WOWA([0.146, 0, 0, 0.5740969929676947, 0, 0.319, 0, 0.696]) # 0.33: 0.365
    print "WOWA test 8: ", fmt.WOWA([1, 0, 0, 0, 1, 0.454, 0, 1]) # 0.971: 1
    print "WOWA test 9: ", fmt.WOWA([0.696, 0, 0, 0, 1, 0.594, 0, 1]) #
    print "WOWA test 10: ", fmt.WOWA([0, 0, 0, 0, 0, 0.006, 0, 1]) # 0.06: 0.199



def test():
    # let's create a fuzzy measure tool to handle e.g. 4 criteria
    fmt1 = FMtool(4, 'Choquet')
    # let's define a new fuzzy measure v by providing directly the values of the subsets of v
    v = [0., 0.15, 0.25, 0.25, 0.35, 0.6, 0.7, 0.75, 0.25, 0.5, 0.6, 0.65, 0.7, 0.85, 0.95, 1.]
    fmt1.setFuzzyMeasure(v)

    # ok, now we need some input data (ie, a vector of criteria)
    z = [0.9, 0.1, 1, 0]
    # in this example, we have 2 criteria that are satisfied (1 and 3)
    # then, we can for instance compute the Choquet integral of this vector using the fuzzy measure defined here above
    print "Choquet integral of vector " + str(numpy.asarray(z)) + " is equal to " + str(fmt1.Choquet(z))

    # now let's try to compute the aggregation of the input vector using an OWA operator
    weights = [.1, .4, .3, .2]
    fmt1.setOWAweights(weights)
    # this weighting vector tends to model a decision-making scheme in which at least 2-3 criteria must be satisfied (no matter which ones),
    # since the first element of the weighting vector has a very low value (0.1)
    print "OWA aggregation using weighting vector " + str(numpy.asarray(weights)) + " is equal to " + str(fmt1.OWA(z))

    # we can also calculate the aggregation of vector 'x' using a generic model, such as "most of..." or "at least k..."
    # the advantage here is that we don't need to specify any weighting vector, as it will be automatically calculated for us
    # the disadvantage is that we have less control over the decision-making model
    print "Aggregation using the model 'most of': " + str(round(fmt1.mostOf(z)))
    print "Aggregation using the model 'at least 2 out of 4':" + str(round(fmt1.atLeast(z, 2)))

    return fmt1


def test2():
    # in this test we compute the Choquet integral by providing the importances (I_i) and interactions indices (I_ij) instead of the the fuzzy measure v
    # we compute then the Choquet integral using formula (2.6) given by Grabsich in [4].
    fmt1 = FMtool(3, 'ChoquetInt')

    x = [.85, .80, 0.10]
    y = [.3, .5, 0.8]

    # let's set the importance factors (Shapley values)
    fmt1.setImportance(1, 0.4)
    fmt1.setImportance(2, 0.4)
    fmt1.setImportance(3, 0.2)

    # there's a synergy of 40% for the criteria pair (1,2)
    fmt1.setInteraction((1, 2), 0.40, 'syn')
    print "Choquet integral of x=%s: %.2f" % (str(x), fmt1.Aggregate(x))
    print "Choquet integral of y=%s: %.2f" % (str(y), fmt1.Aggregate(y))

    # let's aggregate now with OWA and w = [0.1, 0.8, 0.1]
    w = [0.1, 0.8, 0.1]
    fmt1.setOWAweights(w)
    fmt1.setAggrType('OWA')
    print "Using now OWA with w=%s" % str(w)
    print "Aggregated value for x=%s: %.2f" % (str(x), fmt1.Aggregate(x))
    print "Aggregated value for y=%s: %.2f" % (str(y), fmt1.Aggregate(y))

    return fmt1


if __name__ == '__main__':
    test_WOWA()
