from pylab import *
from numpy import random

class RSGenerator(object):

	"""docstring for RSGenerator"""
	def __init__(self, species=1, reactions=1, iterations=1):
		super(RSGenerator, self).__init__()
		self.SPECIES=species
		self.REACTIONS=reactions
		self.ITERATIONS=iterations
		
	def create(self, prefix, binomial_coeff=0.05):

		list_species = map(lambda x: "x"+str(x), xrange(self.SPECIES))
		#print "Species:", list_species

		suffix = "_"+str(self.SPECIES)+"x"+str(self.REACTIONS)+"_bc"+str(binomial_coeff)

		# create the reactions
		reactions=[]
		for reac in xrange(self.REACTIONS):
			n_reac =  binomial(self.SPECIES, binomial_coeff)+1
			n_inhi =  binomial(self.SPECIES, binomial_coeff)
			n_prod =  binomial(self.SPECIES, binomial_coeff)+1
			reactants = random.choice(list_species, size=n_reac, replace=False)
			try:
				inhibitors = random.choice([ x for x in list_species if x not in reactants], size=n_inhi, replace=False)
			except:
				inhibitors = []
			products = random.choice(list_species, size=n_prod, replace=False)
			reactions.append( (reactants, inhibitors, products) )

		with open(prefix+"_gross"+suffix, "w") as fo:
			fo.write("ReactionSystem [ \n\n")
			for (rea, inhi, prod) in reactions:
				re  = "["+",".join(map(lambda x: "\""+x+"\"", rea))+"] "
				ini = " ["+",".join(map(lambda x: "\""+x+"\"", inhi))+"] "
				pro = "["+",".join(map(lambda x: "\""+x+"\"", prod))+"]"
				fo.write("  Reaction "+re+ini+pro+"\n")
			fo.write("\n]")

		with open(prefix+"_orsa"+suffix+".orsa", "w") as fofull:
			with open(prefix+"_brsim"+suffix, "w") as fo:
				for (rea, inhi, pro) in reactions:
					#print rea
					fo.write(" ".join(rea)+", ")
					fo.write(" ".join(inhi)+", ")
					fo.write(" ".join(pro)+"\n")		

					fofull.write(" ".join(rea)+", ")
					fofull.write(" ".join(inhi)+", ")
					fofull.write(" ".join(pro)+"\n")		
			
			fofull.write("\n\n---\n\n")

			# create the context
			with open(prefix+"_gross_context"+suffix, "w") as fo:
				with open(prefix+"_brsim_context"+suffix, "w") as fo2:		
					for it in xrange(self.ITERATIONS):
						n_spec =  binomial(self.SPECIES, 0.4)
						if n_spec==0: 
							fo.write(".\n")
							fo2.write(".\n")
							fofull.write(".\n")
						else:
							species = random.choice(list_species, n_spec, replace=False)
							str_species = map(lambda x: "\""+x+"\"", species)
							fo.write("\t".join(str_species)+"\n")
							fo2.write(" ".join(species)+"\n")
							fofull.write(" ".join(species)+"\n")

NRUNS = 30
NITERS = 1000

if __name__ == '__main__':
        for n in [10,100,1000,2000]:
		for alpha in [0.01,0.05,0.10]:
                        for i in range(NRUNS):
	                        rsg = RSGenerator(species=n, reactions=n, iterations=NITERS)
	                        rsg.create("test_%d" % i, binomial_coeff=alpha)
