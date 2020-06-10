import os
import configparser


class NmrConfig:

    def __init__(self):
        self.autoPlot = True
        self.keepZoom = True
        self.fontSize = 13
        self.phaseReferenceColour = 'Red'
        self.fName = '.pyMetaboLab.config'
        self.homeDir = os.path.expanduser('~')
        self.configFile = os.path.join(self.homeDir, self.fName)
        self.value = ""

    def makeConfig(self):
        config = configparser.ConfigParser()
        autoPlot = 'yes' if self.autoPlot is True else 'no'
        keepZoom = 'yes' if self.keepZoom is True else 'no'
        config['GUI'] = {'autoPlot': autoPlot,
                         'keepZoom': keepZoom,
                         'fontSize': str(self.fontSize)}
        config['Disp'] = {'phaseReferenceColour': self.phaseReferenceColour}
        return config

    def saveConfig(self):
        config = self.makeConfig()
        with open(self.configFile, 'w') as configfile:
            config.write(configfile)

    def readConfig(self):
        config = configparser.ConfigParser()
        config.read(self.configFile)
        if len(config.sections()) == 0:
            config = self.makeConfig()
            self.saveConfig()

        for k in config.sections():
            for l in config[k]:
                self.setValues(l, config[k][l])


    def set_autoplot(self, value):
        self.autoPlot = True if value == "yes" else False

    def set_keepzoom(self, value):
        self.keepZoom = True if value == "yes" else False

    def set_fontsize(self, value):
        self.fontSize = float(value)

    def set_phasereferencecolour(self, value):
        self.phaseReferenceColour = value

    def setValues(self, key, value):
        mName = "self.set_" + key
        eval(mName)(value)

