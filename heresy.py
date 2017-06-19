from PyQt4 import QtGui, QtCore, uic
import sys
import pyRSSIM
import uuid
import os
from subprocess import check_output
import ConfigParser
import time
import iconrc

VERSION = "1.1.0"




class Reaction(object):

	def __init__(self):
		self.name = None
		self.structure = []

	def __repr__(self):
		return self.name +": "+ str(self.structure)




class About(QtGui.QDialog):

	def __init__(self, main_ref=None):
		super(About, self).__init__(None, QtCore.Qt.Window | QtCore.Qt.WindowStaysOnTopHint)
		#self.setWindowModality(QtCore.Qt.WindowModal)
		uic.loadUi('about2.ui', self)
		self.main_ref=main_ref
		#import iconrc
		


class Options(QtGui.QDialog):

	def __init__(self, main_ref=None):
		super(Options, self).__init__(None, QtCore.Qt.WindowStaysOnTopHint)
		uic.loadUi('settings.ui', self)
		self.main_ref=main_ref

	def browse_gpu_file(self):
		fname = QtGui.QFileDialog.getOpenFileName(self, 'Select executable file', '.')
		self.gpusimpath.setText(fname)


	def save_config(self):
		cp = ConfigParser.RawConfigParser()
		cp.add_section('General') 
		self.main_ref.GPUsimulator = str(self.gpusimpath.text())
		cp.set('General', 'GPUsimulator', self.main_ref.GPUsimulator)
		cp.set('General', 'threshold', self.main_ref.THRESHOLD)
		with open("heresyconfig.ini", "wb") as configfile:
			cp.write(configfile)
		print " * Configuration saved in heresyconfig.ini"		
		self.hide()


class MyWindow(QtGui.QMainWindow):


	def __init__(self):
		super(MyWindow, self).__init__()
		uic.loadUi('maingui.ui', self)
		self.show()

		self.heresyfile = None

		self.allowed_species = []

		# defauts 
		self.GPUsimulator = None
		self.THRESHOLD = 128

		# counters
		self.NUMREACTIONS = 0
		self.CONTEXTLENGTH = 0

		# create area for comm
		self.comm = QtGui.QLabel("Retrieving GPU simulator...")
		self.statusBar().addWidget(self.comm,1)
		
		self.entity_comm = QtGui.QLabel("Entities: -")
		self.statusBar().addWidget(self.entity_comm,1)

		self.reactions_comm = QtGui.QLabel("Reactions: -")
		self.statusBar().addWidget(self.reactions_comm,1)

		self.iterations_comm = QtGui.QLabel("Iterations: -")
		self.statusBar().addWidget(self.iterations_comm,1)

		self.flag_normalize.setVisible(False)

		self.setWindowTitle("HERESY v"+VERSION)


	def load_settings(self, option_ref=None):

		# config
		self.load_settings_from_file()

		# update config
		option_ref.threshold.setText(str(self.THRESHOLD))
		
		if self.check_gpu_sim(self.GPUsimulator):			
			self.comm.setText("GPU-powered simulator AVAILABLE")			
			option_ref.gpusimpath.setText(self.GPUsimulator)
		else:
			self.statusBar().showMessage("GPU-powered simulator NOT AVAILABLE")


	def load_settings_from_file(self):
		cp = ConfigParser.RawConfigParser()
		cp.read("heresyconfig.ini")

		try:
			print " * Loading setting for GPU simulator:",
			self.GPUsimulator = cp.get('General', 'GPUsimulator')
			print self.GPUsimulator			
		except ConfigParser.NoSectionError:
			print "WARNING: heresyconfig.ini seems to be broken or missing"
		except ConfigParser.NoOptionError:
			print "WARNING: heresyconfig.ini seems to be old"

		try:
			print " * Loading setting for threshold:",			
			self.THRESHOLD = cp.getint('General', 'threshold')
			print self.THRESHOLD
		except ConfigParser.NoSectionError:
			print "WARNING: heresyconfig.ini seems to be broken or missing"
		except ConfigParser.NoOptionError:
			print "WARNING: heresyconfig.ini seems to be old"


	def check_gpu_sim(self, path):
		try:
			check = os.path.isfile(path) 
			print "Binary file found for GPU simulator"
			return True
		except:
			print "Binary file not found for GPU simulator"
			return False


	def normalization(self, rs, context):

		reaction_number = 0

		species = set()
		reactions = set()
		list_of_reactions = []

		for linea in rs.split("\n"):
			#print linea
			if "Reaction " in linea:
				reaction_number += 1			

				NR = Reaction()
				NR.name = "R"+str(reaction_number)

				linea = linea.strip()
				linea = linea.split()				
				reactants  = eval(linea[1])
				inhibitors = eval(linea[2])
				products   = eval(linea[3])
				NR.structure = [reactants, inhibitors, products]
				for r in reactants: species.add(r)
				for i in inhibitors: species.add(i)
				for p in products: species.add(p)
				#print " Structure:", NR.structure
				reactions.add("R"+str(reaction_number))

				list_of_reactions.append(NR)

		species = list(species)[:]	
		self.allowed_species = species
		reactions = list(reactions)[:]		
		new_species = species[:]
		new_species.extend(reactions)
		new_species.append("HEART")		

		new_reactions = []
		for reaction in list_of_reactions:
		
			for reactant in reaction.structure[0]:	# reactants list
				new_reactions.append([ ['HEART'], [reactant], [reaction.name]])
		
			for inhibitor in reaction.structure[1]:	# inhibitors list
				new_reactions.append([ [inhibitor], ['SPADES'], [reaction.name]])

			for product in reaction.structure[2]:	# products list
				new_reactions.append([ ['HEART'], [reaction.name], [product]])

		#new_reactions.append( [ ['HEART'], [species[0]], ['HEART']] )
		#new_reactions.append( [ ['HEART'], [reactions[0]], ['HEART']] )		
		new_reactions.append( [ ['HEART'], ['SPADES'], ['HEART']] )		
	
		newdesc = ""
		newdesc+= "ReactionSystem [\n\n"
		for nr in new_reactions:
			newdesc+="  Reaction "
			newdesc+= str(nr[0]).replace("'","\"")
			newdesc+= " "
			newdesc+= str(nr[1]).replace("'","\"")
			newdesc+= " "
			newdesc+= str(nr[2]).replace("'","\"")
			newdesc+= "\n"
		newdesc+= "\n]"
		
		newcontext = []
		for x in context.split("\n"):
			#newcontext.append("\"HEART\" "+x)
			newcontext.append(x)
			newcontext.append(".")
		newcontext[0] = "\"HEART\"\t"+newcontext[0]
		newcontext = "\n".join(newcontext)

		return newdesc, newcontext


	def convert_context(self, text):
		text = str(text)
		carriage=text.split("\n")
		conv_text = ""
		for riga in carriage:
			if riga=="": continue
			if riga[0]=="#": continue
			riga = riga.split()
			stringa_ret = ""
			for el in riga:
				if el==".":
					stringa_ret += ".\n"
				else:
					stringa_ret += '"'+el+'"\t'
			conv_text += stringa_ret[:-1]+"\n"
			conv_text = conv_text.replace("'", '"')
		return conv_text


	def convert_description(self, text):
		text = str(text)
		carriage = text.split("\n")
		conv_text = "ReactionSystem [\n\n"
		for riga in carriage:
			if riga=="": continue
			if riga[0]=="#": continue
			riga = riga.split(",")
			
			reagenti = riga[0].split()
			reagenti = "["+",".join(map(lambda x: '"'+x+'"', reagenti))+"]"
			
			inibitori = riga[1].split()
			inibitori = "["+",".join(map(lambda x: '"'+x+'"', inibitori))+"]"
			
			prodotti = riga[2].split()
			prodotti = "["+",".join(map(lambda x: '"'+x+'"', prodotti))+"]"
			
			conv_text += "  Reaction " + reagenti + " "+ inibitori +" "+ prodotti+"\n"
			conv_text = conv_text.replace("'", '"')
		
		conv_text += "\n]"
		return conv_text

	def start_simulation_auto(self, keep_intermediate_files=False):
		desc = self.description.toPlainText()
		#print desc
		desc = self.convert_context(desc)
		detected_rules = len(desc.split("\n"))

		#heuristic
		if detected_rules<128:
			print " * Detected", detected_rules, "rules: simulating with CPU"
			self.start_simulation(keep_intermediate_files=keep_intermediate_files)
		else:
			print " * Detected", detected_rules, "rules: simulating with GPU"
			self.start_simulation_gpu(keep_intermediate_files=keep_intermediate_files)


	def actual_GPU_simulation(self, u1, u2, keep_intermediate_files=False):
		#print " * Launching GPU simulator..."
		u3 = str(uuid.uuid4())
		args = [
			self.GPUsimulator, 
			"-r", u1,
			"-c", u2, 
			"-o", u3
			]

		if self.flag_normalize.isChecked():
			args.append("-l")

		ret = check_output(args)
		ret = ret.replace("\"", "")
		#print ret
		#exit()
		"""
		full_ret = ""
		with open(u3) as fi:
			for ret in fi:		
				#print ret
				ret = ret.replace("\"", "")
				ret = ret.split()
				full_ret += " ".join(ret)+"\n"
		if not keep_intermediate_files:
			os.remove(u3)
		"""
		return ret


	def read_form_and_convert(self):
		desc =  self.description.toPlainText()
		cont =  self.context.toPlainText()

		desc = self.convert_description(desc)
		cont = self.convert_context(cont)

		if self.flag_normalize.isChecked():
			print " * Normalization in progress"
			desc, cont = self.normalization(desc, cont)
			#print desc
			#print cont

		u1 = str(uuid.uuid4())
		u2 = str(uuid.uuid4())
		try:
			with open(u1, "w") as fo: fo.write(desc)
			with open(u2, "w") as fo: fo.write(cont)
		except:
			print "ERROR: cannot create intermediate files"
			exit(-1)

		return desc, cont, u1, u2

	def denormalize_result(self, full_result):
		if self.flag_normalize.isChecked():
			full_result = full_result.split("\n")[::2]
			full_result_clean = []
			#print self.allowed_species
			#print full_result
			for riga in full_result[:-1]:
				riga = riga.split()
				entities_in_result = set(riga[1:])
				full_result_clean.append( str(int(riga[0])/2) + " " + " ".join(entities_in_result.intersection(self.allowed_species)) )
			full_result = "\n".join(full_result_clean)			
		return full_result



	# simulation on GPU
	def start_simulation_gpu(self, keep_intermediate_files=False):
		
		desc, cont, u1, u2 = self.read_form_and_convert()

		start = time.time()
		full_result = self.actual_GPU_simulation(u1, u2, keep_intermediate_files=keep_intermediate_files)
		end = time.time()
		print " * Elapsed time for GPU simulation:", end-start, "seconds"

		#print full_result

		"""
		if self.flag_normalize.isChecked():
			full_result = full_result.split("\n")[::2]
			full_result = "\n".join(full_result)
		
		if self.flag_normalize.isChecked():
			full_result = full_result.split("\n")[::2]
			full_result_clean = []
			for riga in full_result[:-1]:
				riga = riga.split()
				entities_in_result = set(riga[1:])
				full_result_clean.append( str(int(riga[0])/2) + " " + " ".join(entities_in_result.intersection(self.allowed_species)) )
			full_result = "\n".join(full_result_clean)			
		"""
		full_result = self.denormalize_result(full_result)

		self.result.setText( full_result )

		if not keep_intermediate_files:
			os.remove(u1)
			os.remove(u2)


	# simulation on CPU
	def start_simulation(self, keep_intermediate_files=False):
		
		desc, cont, u1, u2 = self.read_form_and_convert()

		RS = pyRSSIM.ReactionSystem()	
		RS.open_system(u1)
		RS.open_context(u2)

		start = time.time()
		full_result = ""
		for n, res in enumerate(RS.simulate_yield()):
			full_result+=str(n)+" "+" ".join(map(str, res))+"\n"

		#print full_result

		"""
		if self.flag_normalize.isChecked():
			full_result = full_result.split("\n")[::2]
			full_result_clean = []
			for riga in full_result[:-1]:
				riga = riga.split()
				entities_in_result = set(riga[1:])
				full_result_clean.append( str(int(riga[0])/2) + " " + " ".join(entities_in_result.intersection(self.allowed_species)) )
			full_result = "\n".join(full_result_clean)			
		"""
		full_result = self.denormalize_result(full_result)

		self.result.setText( full_result )
		end = time.time()
		print " * Elapsed time for CPU simulation:", end-start, "seconds"

		if not keep_intermediate_files:
			os.remove(u1)
			os.remove(u2)


	def validate_description(self, text, verbose=False):
		entities = set()
		text = str(text)
		carriage=text.split("\n")
		numreactions = len(carriage)
		if verbose: print " * Detected", numreactions, "lines"

		carriage = filter(lambda x: x!="", carriage)
		carriage = filter(lambda x: x[0]!="#", carriage)

		numreactions = len(carriage)
		if verbose: print " * Detected", numreactions, "reactions after comments and empty lines removal"
		self.NUMREACTIONS = numreactions

		for riga in carriage:
			riga = riga.split(",")
			if len(riga)!=3: 
				if verbose:
					if len(riga)<3:	print "ERROR: missing comma in reaction", riga
					if len(riga)>3:	print "ERROR: too many commas in reaction", riga
				return False
			else:
				entities.update(riga[0].split())
				entities.update(riga[1].split())
				entities.update(riga[2].split())
			if verbose:
				print len(entities), "entities detected:", entities
			self.entity_comm.setText("Entities: "+str(len(entities)))

		if verbose: print " * All reactions are okay!"
		self.reactions_comm.setText("Reactions: "+str(self.NUMREACTIONS))
		return True

	def validate_context(self, text, verbose=False):
		text = str(text)
		carriage=text.split("\n")
		numreactions = len(carriage)
		if verbose: print " * Detected", numreactions, "lines"

		carriage = filter(lambda x: x!="", carriage)
		carriage = filter(lambda x: x[0]!="#", carriage)

		numsteps = len(carriage)		
		if verbose: print " * Detected", numsteps, "contextes after comments and empty lines removal"
		self.CONTEXTLENGTH = numsteps

		if numsteps==0:
			#print "ERROR: no context detected"
			return False
		
		if verbose: print " * Context seems okay!"
		self.iterations_comm.setText("Iterations: "+str(self.CONTEXTLENGTH))
		return True


	def try_unlock(self):
		desc =  self.description.toPlainText()
		cont =  self.context.toPlainText()

		pal = QtGui.QPalette()
		okay = QtGui.QColor(200, 200, 255)
		wrong = QtGui.QColor(255, 200, 200)

		descok=self.validate_description(desc)
		contok=self.validate_context(cont)
		
		if descok:
			pal.setColor(QtGui.QPalette.Base, okay)			
		else:
			pal.setColor(QtGui.QPalette.Base, wrong)
		self.description.setPalette(pal)

		if contok:
			pal.setColor(QtGui.QPalette.Base, okay)			
		else:
			pal.setColor(QtGui.QPalette.Base, wrong)		
		self.context.setPalette(pal)

		self.button_simulate.setEnabled(contok and descok)
		self.simwithgpu.setEnabled(contok and descok)
		self.automatic.setEnabled(contok and descok)


	def open_preferences(self):
		options_window.show()

	def updated_reactions(self):
		t  = self.description.toPlainText()
		#print t
		#print " * Reactions valid:", self.validate_description(t), self.NUMREACTIONS

	def updated_contexts(self):
		t  = self.context.toPlainText()
		#print " * Context valid:", self.validate_context(t), self.CONTEXTLENGTH


	def show_about(self):
		about_window.show()


	def saveas(self):
		fname = QtGui.QFileDialog.getSaveFileName(self, 'Select target HERESY file', '.', '*.rsy')
		if fname!="":
			self.save_heresy_file(fname)

	def save(self):
		if self.heresyfile == None:
			self.saveas()
		else:
			self.save_heresy_file(self.heresyfile)

	def save_heresy_file(self, targetfile):
		with open(targetfile, "w") as fo:
			fo.write(self.description.toPlainText())
			fo.write("\n\n---\n\n")
			fo.write(self.context.toPlainText())
		print " * File saved to", targetfile
		self.heresyfile = targetfile
		self.setWindowTitle("HERESY v"+VERSION+" - "+fname)

	def openfile(self):
		fname = QtGui.QFileDialog.getOpenFileName(self, 'Select HERESY file', '.', '*.rsy')
		if fname!="":
			self.open_heresy_file(fname)
			self.setWindowTitle("HERESY v"+VERSION+" - "+fname)


	def open_heresy_file(self, fname):
		with open(fname) as fi:
			A = fi.readlines()
			index = 0
			for a in A:				
				if a[:3] == "---": break
				index+=1				
			self.description.setPlainText("".join(A[0:index-2]))
			self.context.setPlainText("".join(A[index+2:]))
			self.heresyfile = fname
			print " * HERESY file loaded:", fname



if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)
    window = MyWindow()
    window.show()
    options_window = Options(main_ref=window)    
    about_window = About(main_ref=window)

    window.load_settings(option_ref=options_window) 

    sys.exit(app.exec_())
