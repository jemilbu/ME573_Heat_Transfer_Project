import csv
import matplotlib.pyplot as plt
import matplotlib as mpl

class View():
    def __init__(self):
        print("  Welcome to the Project for Team 1  ".center(80,'#'))
        print("  Team Members:  ".center(80,'#'))
        print("  Jacob Milburn, Nathan Lauritsen, Isaiah Elsner,  ".center(80))
        print("  Braden Harter and Daniel Lopez  \n".center(80))
        self.fileName()

    def fileName(self):
        while True:
            try:
                inCSV = input("Do you want to save as a CSV (y/n): ")
                if (inCSV[0] == 'y'): self.runCSV = True
                elif (inCSV[0] == 'n'): self.runCSV = False
                else:
                    print("Invalid input, please try again\n")
                    continue
            except:
                print("Invalid input, please try again\n")
                continue
            break

        while True:
            try:
                inPlot = input("Do you want to save as a plot (y/n): ")
                if (inPlot[0] == 'y'): self.runPlot = True
                elif (inPlot[0] == 'n'): self.runPlot = False
                else:
                    print("Invalid input, please try again\n")
                    continue
            except:
                print("Invalid input, please try again\n")
                continue
            break
        
        if (self.runCSV == True or self.runPlot == True):
            self.filename = input("Please enter an output file name: ")

    def askToSee(self,models,mdot):
        while True:
            try:
                inSee = input("Do you want to see a specific mass flow rate? (y/n): ")
                if (inSee[0] == 'y'): self.runSee = True
                elif (inSee[0] == 'n'): self.runSee = False
                else:
                    print("Invalid input, please try again\n")
                    continue
            except:
                print("Invalid input, please try again\n")
                continue
            
            if (self.runSee == True):
                while True:
                    while True:
                        try:
                            inSee = input("Enter a mass flow rate in kg/hr: ")
                            if ',' in inSee:
                                print("Please enter the value again with no commas\n")
                                continue
                            self.mdotCall = int(inSee)
                            for i in range(len(mdot)):
                                if mdot[i] == self.mdotCall:
                                    print("")
                                    self.printValues(models[i].giveAll())
                                    foundFlag = True
                                    break
                                else:
                                    foundFlag = False
                            if foundFlag:
                                foundFlag = False
                                break
                            else:
                                print("Couldn't find that value, please try again\n")
                                continue
                        except:
                            print("Invalid input, please try again\n")
                            continue
                        break
                    
                    while True:
                        try:
                            runAgian = input("Do you want to go again? (y/n): ")
                            if (runAgian[0] == 'y'): 
                                self.runAgainFlag = True
                            elif (runAgian[0] == 'n'): 
                                self.runAgainFlag = False
                            else:
                                print("Invalid input, please try again\n")
                                continue
                        except:
                            print("Invalid input, please try again\n")
                            continue
                        break

                    if self.runAgainFlag == True:
                        continue
                    else:
                        break
            break

    def printValues(self, vals):
        print("Mass Flow Rate: {:.2f}".format(vals[0]))
        print("Tcold Out: {:.4f}".format(vals[1]))
        print("Thot Out: {:.4f}".format(vals[2]))
        print("Cc: {:.4f}".format(vals[3]))
        print("Ch: {:.4f}".format(vals[4]))
        print("Re: {:.4f}".format(vals[5]))
        print("ΔP: {:.4f}".format(vals[6]))
        print("Nu: {:.4f}".format(vals[7]))
        print("hc (water): {:.4f}".format(vals[8]))
        print("Thermal Resistance: {:.4f}".format(vals[9]))
        print("U: {:.4f}".format(vals[10]))
        print("NTU: {:.4f}".format(vals[11]))
        print("Effectivness: {:.4f}".format(vals[12]))
        print("Q : {:.4f}".format(vals[13]))
        print("")

    def outputTableToCSV(self,vals):
        if (self.runCSV == True):
            with open(self.filename+".csv",'w', newline='') as csvfile:
                outfile = csv.writer(csvfile, delimiter=',')
                outfile.writerow(["Mass Flow Rate (kg/hr)","T_cold_out (oC)", 
                "T_hot_out (oC)","Cc (W/K)","Ch (W/K)","ReD (Water)",
                "P (Pressure Drop) (N/m2)","NUD (Water)","hc (Water) (W/m^2*K)",
                "Thermal Resistance R (K/W)","Uh (Oil) (W/m^2*K)","NTU",
                "E (Effectiveness)","q (Heat Transfer Rate) (W)"])
                for i in vals:
                    a = i.giveAll()
                    outfile.writerow(a)

    def makePlot(self,mdot,h,tco,tho,e):
        if (self.runPlot == True):
            mpl.rc('font',family='serif',size=10)
            fontTitle = {'family': 'serif','size':'15'}

            markerSize = 3
            plt.plot(mdot,tco,'o--',label="Tc Out (K)",ms=markerSize)
            plt.plot(mdot,tho,'o-.',label="Th Out (K)",ms=markerSize)
            plt.plot(mdot,h,'o:',label="H (W/m^2 K) (1x10^1)",ms=markerSize)
            plt.plot(mdot,e,'o-',label="% Efficiency (1x10^-1)",ms=markerSize)
            plt.legend()
            plt.title("Heat Transfer Project Results for Team 1",fontdict=fontTitle)
            plt.xlabel("Mass Flow Rate (kg/hr)")
            plt.grid(True)
            plt.savefig(self.filename+"Plot.png")
            plt.show()
        

