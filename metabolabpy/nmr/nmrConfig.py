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
        self.posCol10 = 0.0
        self.posCol11 = 0.0
        self.posCol12 = 1.0
        self.negCol10 = 1.0
        self.negCol11 = 0.0
        self.negCol12 = 0.0
        self.posCol20 = 0.8
        self.posCol21 = 0.8
        self.posCol22 = 1.0
        self.negCol20 = 1.0
        self.negCol21 = 0.8
        self.negCol22 = 0.8
        self.mode = 'light'

    def makeConfig(self):
        config = configparser.ConfigParser()
        autoPlot = 'yes' if self.autoPlot is True else 'no'
        keepZoom = 'yes' if self.keepZoom is True else 'no'
        config['GUI'] = {'autoPlot': autoPlot,
                         'keepZoom': keepZoom,
                         'fontSize': str(self.fontSize),
                         'mode': self.mode}
        config['Disp'] = {'phaseReferenceColour': self.phaseReferenceColour,
                          'posCol10': self.posCol10,
                          'posCol11': self.posCol11,
                          'posCol12': self.posCol12,
                          'negCol10': self.negCol10,
                          'negCol11': self.negCol11,
                          'negCol12': self.negCol12,
                          'posCol20': self.posCol20,
                          'posCol21': self.posCol21,
                          'posCol22': self.posCol22,
                          'negCol20': self.negCol20,
                          'negCol21': self.negCol21,
                          'negCol22': self.negCol22}
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

    def set_poscol10(self, value):
        self.posCol10 = float(value)

    def set_poscol11(self, value):
        self.posCol11 = float(value)

    def set_poscol12(self, value):
        self.posCol12 = float(value)

    def set_negcol10(self, value):
        self.negCol10 = float(value)

    def set_negcol11(self, value):
        self.negCol11 = float(value)

    def set_negcol12(self, value):
        self.negCol12 = float(value)

    def set_poscol20(self, value):
        self.posCol20 = float(value)

    def set_poscol21(self, value):
        self.posCol21 = float(value)

    def set_poscol22(self, value):
        self.posCol22 = float(value)

    def set_negcol20(self, value):
        self.negCol20 = float(value)

    def set_negcol21(self, value):
        self.negCol21 = float(value)

    def set_negcol22(self, value):
        self.negCol22 = float(value)

    def set_mode(self, value):
        self.mode = value

    def setValues(self, key, value):
        mName = "self.set_" + key
        eval(mName)(value)

