import pylab

inputDataFile = 'output.dat'
data = pylab.load(inputDataFile, comments='!', skiprows=0)

l = pylab.plot(data[:,1], 'r-', linewidth=1)

pylab.xlabel('Iteration')
pylab.ylabel('Mean Square Metric')
pylab.title('Registration')
pylab.grid(True)
pylab.savefig('outputMetric.png')

pylab.show()

l = pylab.plot(data[:,2], data[:,3], 'r-', linewidth=1)

pylab.xlabel('X')
pylab.ylabel('Y')
pylab.title('Registration')
pylab.grid(True)
pylab.savefig('outputTranslation.png')

pylab.show()