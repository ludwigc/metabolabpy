import os
import configparser


class NmrConfig:
    
    def __init__(self):
        self.autoPlot                 = True
        self.keepZoom                 = True
        self.fontSize                 = 13
        self.phaseReferenceColour     = 'Red'
        self.fName                    = '.pyMetaboLab.config'
        self.homeDir                  = os.path.expanduser('~')
        self.configFile               = os.path.join(self.homeDir, self.fName)
        
    def makeConfig(self):
        config            = configparser.ConfigParser()
        autoPlot          = 'yes'
        if(self.autoPlot==False):
            autoPlot = 'no'
            
        keepZoom = 'yes'
        if(self.keepZoom==False):
            keepZoom = 'no'
        config['GUI']  = {'autoPlot': autoPlot,
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
        if(len(config.sections())==0):
            config = self.makeConfig()
            self.saveConfig()
        
        self.autoPlot = True
        try:
            if(config['GUI']['autoPlot'] == 'no'):
                self.autoPlot = False

        except:
            pass

        self.keepZoom = True
        try:
            if(config['GUI']['keepZoom'] == 'no'):
                self.keepZoom = False
                
        except:
            pass

        try:
            self.fontSize                 = float(config['GUI']['fontSize'])
        except:
            pass
        
        try:
            self.phaseReferenceColour     = config['Disp']['phaseReferenceColour']
        except:
            pass





