from copy import deepcopy


class Species(object):
	"""docstring for Species"""
	def __init__(self, name):
		super(Species, self).__init__()
		self.name=name
	
	def __repr__(self):
		return "<S. "+str(self.name)+">"

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self, other):
		if isinstance(other, Species):
			return self.__hash__() == other.__hash__()
		else:
			return False

class Reaction(object):
	"""docstring for Reaction"""
	def __init__(self, reactants, inhibitors, products, name=None):
		super(Reaction, self).__init__()
		self.reactants  = reactants
		self.inhibitors = inhibitors
		self.products   = products
		if name==None:
			self.reaction_name = "unnamed"
		else:
			self.reaction_name = name		
		
	def __str__(self):
		return "<Reaction "+self.reaction_name+": "+str(self.reactants)+"->"+str(self.products)+" inhibited by "+str(self.inhibitors)+">"

	def __repr__(self):
		return self.reaction_name
		

class ReactionSystem(object):

	"""docstring for ReactionSystem"""
	def __init__(self):
		super(ReactionSystem, self).__init__()
		self.reactions = set()
		self.species   = set()
		self.state 	   = set()
		self.next_state= set()
		self.interleaved = False
		#print " * New empty Reaction System created"
		
	def simulation_step(self, verbose=False):

		#print " * State right now:", self.state
		for reaction in self.reactions:
			if verbose: print "   Applying reaction:", reaction
			if verbose: print "    Verifying reactants:", reaction.reactants
			if set(reaction.reactants).issubset(self.state):
				if verbose: print "     Reactants are present, veryfing inhibitors:", reaction.inhibitors
				if set(reaction.inhibitors).intersection(self.state)!=set([]) and len(reaction.inhibitors)>0:
					if verbose: print "     Inhibitors are present, reaction aborted"
					# if verbose: print set(reaction.inhibitors).intersection(self.state)
				else:
					if verbose: print "     Inhibitors are not present, reaction will fire adding species:", reaction.products
					for el in reaction.products:
						self.next_state.add(el)
			else:
				if verbose: print "     Reactants not found"
			if verbose: print 
		if verbose: print " > Eventually, the new state is:", self.next_state
		#; print "*"*100




	def simulate(self, steps=None, verbose=False):
		if steps==None:
			steps = len(self.context)-1
			if verbose: print " * Simulation will perform", steps, "steps"			
			if verbose: print " * Results sequence:"

		if verbose: print map(lambda x: x.name, self.state)
		for step in xrange(steps):			
			# print " * Simulation step", step
			self.next_state.clear()
			self.simulation_step(verbose=verbose)
			for el in self.context[step+1]:
				self.next_state.add(el)
			self.state = deepcopy(self.next_state)
			if verbose: print map(lambda x: x.name, self.state)


	def simulate_yield(self, steps=None, verbose=False):
		if steps==None:
			steps = len(self.context)
			if verbose:
				print " * Simulation will perform", steps, "steps"			
				print " * Starting configuration:", self.state

		yield map(lambda x: x.name, self.state)
		for step in xrange(steps):			
			if verbose: print " * Simulation step", step
			self.next_state.clear()
			self.simulation_step(verbose=verbose)
			if (not self.interleaved) or (self.interleaved and step%2==0): 
				yield map(lambda x: x.name, self.next_state)
			#if (not self.interleaved) or (self.interleaved and step%2==0): 
			if step<steps-1:
				if verbose:
					print " >> Applying context:", self.context[step+1]
				for el in self.context[step+1]:
					self.next_state.add(el)

				if verbose:
					print " >>> Finally, new state is", self.next_state
					print "*"*100
			self.state = deepcopy(self.next_state)
			

	def set_initial_state(self, state):
		self.state = deepcopy(state)
		# print " * Initial state set to:", self.state

	def add_reaction(self, reaction):
		self.reactions.add(reaction)
		#print reaction, "added"

	# this method opens a RS encoded in our proprietary format (not brsim)
	def open_system(self, path, verbose=False):
		try:
			with open(path) as fi:
				for linea in fi:
					if linea.find("Reaction ")>0:
						#print linea.split()[1:]						

						reacts, inhibs, prods = map(eval, linea.split()[1:])
						if isinstance(prods, tuple):
							prods = prods[0]
						list_reactants  = []
						list_inhibitors = []
						list_products   = []
						for reactant in reacts:
							list_reactants.append(Species(reactant))
							self.species.add(Species(reactant))
						for inhibitor in inhibs:
							list_inhibitors.append(Species(inhibitor))
							self.species.add(Species(inhibitor))
						for product in prods:
							list_products.append(Species(product))
							self.species.add(Species(product))
						# print "Reaction:", list_reactants, list_inhibitors, list_products
						# reacts = eval(reacts)
						# print reacts, inhibs, prods
						r = Reaction( list_reactants, list_inhibitors, list_products, name="R"+str(len(self.reactions)))
						self.add_reaction(r)	
						# if verbose: print r

					else:
						pass
						# print linea
				
		except IOError:
			print "ERROR: cannot open", path
			exit(-1)


	def open_context(self, path, verbose=False):
		self.context = []
		with open(path) as fi:
			for linea in fi:
				if linea[0]==".":
					self.context.append([])
				else:
					self.context.append(map(lambda x: Species(x), map(eval, linea.split())))
		if verbose: print " * Loaded context:", self.context
		self.rewind()

	def rewind(self):
		self.set_initial_state( self.context[0] )		

	"""
	def renormalize(self, verbose=True):
		# Converts a RS into its minimal normal form. 
		
		A_x = set()
		S_x = deepcopy(self.species)
		S_x.add( Species("HEART"))
		for a in self.reactions: 
			S_x.add( Species(a.__repr__()) )

		for a in self.reactions:
			if verbose: print a
			for r in a.reactants:
				if verbose: print r, 
				# new reaction of type #6
				new_r = Reaction( [Species("HEART")], [r], [Species(a.__repr__())], name="R"+str(len(A_x)) )
				A_x.add(new_r)
				if verbose: print new_r
			
			for i in a.inhibitors:
				if verbose: print i, 
				# new reaction of type #3
				new_r = Reaction( [i], [Species(a.__repr__())], [Species(a.__repr__())], name="R"+str(len(A_x)) )
				A_x.add(new_r)
				if verbose: print new_r
			
			for p in a.products:
				if verbose: print p, 
				# new reaction of type #4
				new_r = Reaction( [Species("HEART")], [Species(a.__repr__())], [p], name="R"+str(len(A_x)) )
				A_x.add(new_r)
				if verbose: print new_r
			if verbose: print 
		if verbose: print "RENORMALIZATION COMPLETED"
		if verbose: print "*"*100	
		if verbose: print " S =", self.species
		if verbose: print " A =", self.reactions
		if verbose: print " S'=", S_x
		if verbose: print " A'=", A_x

		t = list(self.species)[0]
		b = list(self.reactions)[0]

		A_x.add ( Reaction([Species('HEART')], [t], [Species('HEART')]) )
		A_x.add ( Reaction([Species('HEART')], [b.__repr__()], [Species('HEART')]) )

		self.species = S_x
		self.reactions = A_x
	"""

	def set_interleaved_execution(self, interleaved):
		self.interleaved = interleaved
		if interleaved:
			print " * Interleaved execution enabled"


if __name__ == '__main__':
	
	input_path = "C:\\Users\\aresio\\Documents\\Coding\\SimReacSystem\\SimReacSystem\\prova_gross"
	context_path = "C:\\Users\\aresio\\Documents\\Coding\\SimReacSystem\\SimReacSystem\\prova_gross_context"

	RS = ReactionSystem()
	RS.open_system(input_path, verbose=True)
	RS.open_context(context_path, verbose=True)

	for res in RS.simulate_yield(verbose=False):
		print res,
	print 

	RS.renormalize(verbose=False)
	RS.rewind()
	RS.set_interleaved_execution(False)

	for res in RS.simulate_yield(verbose=False):
		print res,

