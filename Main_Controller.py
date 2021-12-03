from Program_Model import Model
from User_View import View

def main():
    v = View()

    m = []
    mdot = []
    h = []
    tco = []
    tho = []
    e = []
    for i in range(5,16):
        m.append(Model(i*1000))
        j = i - 5
        mdot.append(m[j].giveMdot())
        h.append(m[j].giveH() / 10)
        tco.append((m[j].giveTcout()))
        tho.append((m[j].giveThout()))
        e.append(m[j].giveEff() * 1000)

    
    v.makePlot(mdot, h, tco, tho, e)
    v.outputTableToCSV(m)
    v.askToSee(m,mdot)

if __name__ == '__main__':
    main()