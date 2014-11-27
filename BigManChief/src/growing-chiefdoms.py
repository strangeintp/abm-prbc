'''
Created on Nov 22, 2014

@author: stran_000
'''

from random import random
from random import randint
from random import choice
from random import shuffle
from random import uniform
from math import exp
from math import log
from math import sin
from math import pi
from utility import *
from copy import deepcopy

PERIOD = 30
MAX_EFF_TEMP = 25
BIOMASS_SCALE = 1000*2000*PERIOD #pop density 1000 * kCal/day * days/period / max effective temperature
OPTIMAL_TASK_GROUP_SIZE = 6  #scalar stress
MAX_KIN_SELECTION = 4
MAX_AGE = 80
AGE_ADULTHOOD = 20
AGE_SENIOR = 50
MAX_HEALTH = 1.25
WORLD_DIM = 50
TheWorld = None
ThePeople = None
INITIAL_POPULATION = 60
BIRTH_RATE = 0.25 # = 1 - (1-rate_month)**12
FERT_PROB = 1 - (1-BIRTH_RATE)**(1/12)
CONSTRAINED_WORLD = False

class Ancestry(object):
    
    def __init__(self, person, generations=2):
        self.personID = person.ID #the actual person object represented by this node
        self.ancestors = [None, None]
        self.past_length = generations  # how far back the person tracks ancestors
        parents = person.parents
        if len(parents)==2:
            self.ancestors[FEMALE] = [parents[FEMALE].ID, Ancestry.Generate(parents[FEMALE].ancestry.ancestors, generations-1)]
            self.ancestors[MALE] = [parents[MALE].ID, Ancestry.Generate(parents[MALE].ancestry.ancestors, generations-1)]
    
    @staticmethod
    def Generate(parent_ancestors, generations):
        #recursively add parent ancestries until generation tracking variable equals 0
        ancestors = [None, None]
        if generations>0:
            if ancestors[0] != None:
                ancestors[0] = [parent_ancestors[0][0], Ancestry.Generate(parent_ancestors[0][1], generations-1)]
            if ancestors[1] != None:
                ancestors[1] = [parent_ancestors[1][0], Ancestry.Generate(parent_ancestors[1][1], generations-1)]
        return ancestors
    
    def getAncestors(self, ancestors=None, generations=-1):
        if generations==-1:
            generations = self.past_length
        if ancestors==None:
            ancestors = self.ancestors
        ancestor_set = set()
        if generations<=0:
            return ancestor_set
        if ancestors[0] != None:
            ancestor_set = ancestor_set.union(set([ancestors[0][0]]))
            ancestor_set = ancestor_set.union(self.getAncestors(ancestors[0][1], generations-1))
        if ancestors[1] != None:
            ancestor_set = ancestor_set.union(set([ancestors[1][0]]))
            ancestor_set = ancestor_set.union(self.getAncestors(ancestors[1][1], generations-1))
        return ancestor_set
              
    def isKin(self, other):
        my_ancestors = self.getAncestors() | set([self.personID])
        their_ancestors = other.getAncestors() | set([other.personID])
        return my_ancestors & their_ancestors


class Patch(object):
    
    def __init__(self, x, y, avg_annual_biomass, peak_ampl, phase):
        """
        x,y = patch location
        avg_annual_biomass = average annual biomass of the patch; should be a function of effective temperature
        initial_amplitude = multiplies avg_annual_biomass for starting value
        """ 
        self.x = x
        self.y = y
        self.avg_annual_biomass = avg_annual_biomass
        self.peak = peak_ampl
        self.phase = phase
        amplitude = 1.0 + self.peak*sin(-self.phase)
        self.current_biomass = amplitude*avg_annual_biomass
    
    def setAverage(self, average_biomass):
        self.avg_annual_biomass = average_biomass
        
    def setPeak(self, peak_ampl):
        self.peak = peak_ampl
        
    def setPhase(self, phase):
        self.phase = phase
        
    def deplete(self, amount):
        self.current_biomass -= amount
        
    def replenish(self, t):
        """
        """
        prev_ampl = 1.0 + self.peak*sin(2*pi*((t-1)/12)-self.phase)
        amplitude = 1.0 + self.peak*sin(2*pi*(t/12)-self.phase)
        """
        The following rate calculation forces the logistic equation to follow the
        sinusoidal variation in amplitude in the absence of resource depletion
        (e.g., as a result of man's foraging activities)
        """
        rate = amplitude/prev_ampl
        capacity = amplitude*self.avg_annual_biomass
        self.current_biomass += rate*self.current_biomass*(1 - self.current_biomass/capacity)
        
    def forage(self, expenditure):
        r = self.current_biomass/BIOMASS_SCALE
        avg_return = GenBoundedRandomNormal(r, r/3, 0, 2*r)
        amount = expenditure*avg_return
        if amount > self.current_biomass:
            amount = 0.9 * self.current_biomass
        self.deplete(amount)
        return amount
    
    def __sub__(self, other):
        dx = self.x-other.x
        dy = self.y-other.y
        if not CONSTRAINED_WORLD:
            dx %= int(WORLD_DIM/2)
            dy %= int(WORLD_DIM/2)
        return (dx**2 + dy**2)**0.5
            
        
class World(object):
    
    @staticmethod
    def generateWithHomogenousPatches(world, dimension, eff_temp):
        patch_biomass = BIOMASS_SCALE*eff_temp/MAX_EFF_TEMP
        peak_ampl = (MAX_EFF_TEMP-eff_temp)/MAX_EFF_TEMP
        for x in range(dimension):
            for y in range(dimension):
                patch = Patch(x, y, patch_biomass, peak_ampl, 0.0)
                world.patches[(x,y)] = patch
                world.people_at[patch] = []
    
    @staticmethod
    def generate(world, dimension, eff_temp, spatial_variance=0, temporal_variance=0):
        World.generateWithHomogenousPatches(world, dimension, eff_temp)
        x_vals = [i for i in range(WORLD_DIM)]
        y_vals = [i for i in range(WORLD_DIM)]
        shuffle(x_vals)
        shuffle(y_vals)
        locations = list(zip(x_vals,y_vals))
        #increase half the cells
        for (x,y) in locations[0:int(WORLD_DIM/2)]: 
            here = world.patches[(x,y)]
            neighbors = world.neighborsOf(here)
            neighbors.remove(here)
            inc = here.avg_annual_biomass*spatial_variance
            here.setAverage(here.avg_annual_biomass + inc)
            if temporal_variance > 0:
                phase = int(GenBoundedRandomNormal(0, temporal_variance, -4.5, 4.5))
                here.setPhase(phase)
            for there in neighbors: #spread the increment around
                dec = there.avg_annual_biomass - inc/8
                there.setAverage(dec)
        #decrease the other half
        for (x,y) in locations[int(WORLD_DIM/2):]: 
            here = world.patches[(x,y)]
            neighbors = world.neighborsOf(here)
            neighbors.remove(here)
            dec = here.avg_annual_biomass*spatial_variance
            here.setAverage(here.avg_annual_biomass - dec)
            if temporal_variance > 0:
                phase = int(GenBoundedRandomNormal(0, temporal_variance, -4.5, 4.5))
                here.setPhase(phase)
            for there in neighbors: #spread the decrement around
                inc = there.avg_annual_biomass + dec/8
                there.setAverage(inc)
        
    
    def __init__(self, dimension, generator, eff_temp=25, spatial_variance=0, temporal_variance=0):
        print("Generating the world...")
        self.patches = {}
        self.people_at = {}
        self.eff_temp = eff_temp
        generator(self, dimension, eff_temp=eff_temp, spatial_variance=spatial_variance, temporal_variance=temporal_variance)
        
    def placePersonAt(self, person, patch):
        self.people_at[patch].append(person)
        person.location = patch
        
    def removePersonFrom(self, person, patch):
        self.people_at[patch].remove(person)
        
    def movePersonTo(self, person, patch):
        self.removePersonFrom(person, person.location) 
        self.placePersonAt(person, patch)
        
    def step(self, t):
        for patch in self.patches.values():
            patch.replenish(t)
        self.month = t%12
            
    def neighborsOf(self, patch):
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                x = patch.x+dx
                y = patch.y+dy
                if not CONSTRAINED_WORLD:
                    neighbors.append(self.patches[(x%WORLD_DIM,y%WORLD_DIM)])
                elif x>0 and y<0 and x<WORLD_DIM and y<WORLD_DIM:
                    neighbors.append(self.patches[(x,y)])
        shuffle(neighbors)
        return neighbors
        

"""
Indeces for attributes
"""
coll = 0 # what fraction of energy does an agent put into collective action
prestige_weight = 1 # for annealing; how much does an agent value prestige
corp_prestige = 2 # prestige gained corporate or network; values [0,1]
kin_selection = 3
grp2size = 4
cultural_attributes = [coll, prestige_weight, corp_prestige, kin_selection, grp2size]
attribute_ranges = [1, 1, 1, MAX_KIN_SELECTION, 500]
MALE = 1
FEMALE = 0

BOILING = 100
delta_E = (len(cultural_attributes)*0.25)**0.5
k = -delta_E/(BOILING*log(0.5))
mu = 0.01
def anneal(group):
    temperature = 0
    current_energy = 0
    avg_coll = 0
    for person in group:
        temperature += person.getTemperature(group) + 1/len(group)
        avg_coll += person.attributes[coll]
        current_energy += person.calculateEnergy(group)
    if not group:
        print("no people!")
    temperature /= len(group)
    current_energy /= 2  # corrects for double-counting in loop above
    avg_coll /= len(group)
    #print("T = %f \t E = %f \t coll = %f"%(temperature, current_energy, avg_coll))
    #print("patch resource = %f"%(group[0].location.current_biomass/BIOMASS_SCALE))
    group = sorted(group, key = lambda person : person.temperature, reverse = True)
    for person in group:  # start annealing from the highest temperature agent
        new_person = person.mutate(group)
        person.calculateEnergy(group)
        delta_energy = new_person.energy - person.energy
        if delta_energy <= 0:  # accept the new person if lower energy
            person.attributes = new_person.attributes
        elif random() < exp(-delta_energy/(k*temperature)): # or accept with some probability
            person.attributes = new_person.attributes
        # otherwise, discard the "new" person
    return temperature
        
def synergize(n):
    load = (n**2 - n)/2
    n_best = OPTIMAL_TASK_GROUP_SIZE
    capacity = (n_best**2 - n_best)/2
    if n<n_best:
        synergy = 1 + load/capacity
    else:
        synergy = 1 + capacity/load
    return synergy

class Person(object):
    
    nextID = -1
    @staticmethod
    def getNextID():
        Person.nextID += 1
        return Person.nextID
    
    @staticmethod
    def generateRandomAdult():
        person = Person()
        person.attributes[coll] = random()
        person.attributes[prestige_weight] = random()
        person.attributes[corp_prestige] = random()
        person.attributes[kin_selection] = randint(1, MAX_KIN_SELECTION)
        person.attributes[grp2size] = randint(25, 75)
        person.age -= 1/12
        person.grow()
        person.ancestry = Ancestry(person, person.attributes[kin_selection])
        return person
    
    @staticmethod
    def generateFromParents(mother, father):
        person = Person(parents=[mother,father])
        person.age = -1/12
        person.grow()
        if person.gender==MALE:
            parent0 = father
            parent1 = mother
        else:
            parent0 = mother
            parent1 = father
        for attribute in cultural_attributes:
            beta = GenBoundedRandomNormal(0.25, 0.1, -0.25, 0.5)
            value = (1-beta)*parent0.attributes[attribute] + (beta)*parent1.attributes[attribute]
            value = max(0, value)
            value = min(value, attribute_ranges[attribute])
            person.attributes[attribute] = value
        
        # foraging ability
        beta = GenBoundedRandomNormal(0.25, 0.1, -0.25, 0.5)
        person.ability = (1-beta)*parent0.ability + (beta)*parent1.ability
        
        person.location = mother.location  
        person.ancestry = Ancestry(person, person.attributes[kin_selection])      
        
        return person
    
    def __init__(self, parents = []):
        self.ID = Person.getNextID()  
        self.attributes = []
        self.age = randint(AGE_ADULTHOOD, AGE_SENIOR)
        self.health = 1.0
        self.gender = MALE if random()>0.5 else FEMALE
        self.ability = GenBoundedRandomNormal(1.0, 0.2, 0.75, 1.25)
        self.current_ability = self.ability
        self.temperature = 0
        self.location = None
        self.food_stored = 0
        for attribute in cultural_attributes:
            self.attributes.append(0)
        self.spouses = []
        self.children = []
        self.parents = parents
        self.group1 = None
        self.group2 = None
        self.group3 = None
        self.pregnancy = 0
        self.caloric_need = 2000*30
        self.family = Family()
        self.patch_memory = None
        
    def dies(self):
        for spouse in self.spouses:
            spouse.spouses.remove(self)
        for child in self.children:
            if self in child.parents:
                child.parents.remove(self)
        for parent in self.parents:
            if self in parent.children:
                parent.children.remove(self)
    
    def getTemperature(self, group):
        self.temperature = self.getHealthTemperature()        
        # self.temperature += self.getReproductiveTemperature()
        return self.temperature
    
    def getHealthTemperature(self):
        return BOILING*4*(MAX_HEALTH - self.health)/3.0  # scales to 0 at health = 1.0 and 100 at health=0.25
    
    def getPrestigeTemperature(self):
        """
        returns an agent's social temperature based on prestige among peer group
        """
        return 0 # stub value
    
    def culturalDistanceFrom(self, other):
        distance = sum([((self.attributes[attribute] - other.attributes[attribute])/attribute_ranges[attribute])**2 for attribute in cultural_attributes])
        return distance**(0.5)
    
    def calculateEnergy(self, group):
        self.energy = 0
        for person in group:
            if not person==self:
                self.energy += self.culturalDistanceFrom(person)
            
        return self.energy      
    
    def mutate(self, group):
        """
        mutates an individual's cultural attributes
        moves attributes randomly towards other random members
        """
        new_person = deepcopy(self)
        for attribute in cultural_attributes:
            other = choice(group)
            mu_local = GenBoundedRandomNormal(0, mu, -3*mu, 3*mu)
            new_attr = (1-mu)*self.attributes[attribute] + mu*other.attributes[attribute]
            new_attr = max(0, new_attr)  # clamp new attribute value to valid range
            new_attr = min(new_attr, attribute_ranges[attribute])
            new_person.attributes[attribute] = new_attr
        new_person.calculateEnergy(group)
        
        return new_person
    
    def shareWithNuclearFamily(self, amount):
        family = [self]+self.spouses+self.children
        capacities = [(MAX_HEALTH-member.health)*member.caloric_need for member in family]
        total_cap = sum(capacities)
        portions = [amount*capacity/total_cap for capacity in capacities]
        for member, portion in zip(family, portions):
            member.eat(portion)
    
    def eat(self, amount):
        food_capacity = (MAX_HEALTH-self.health)*self.caloric_need
        if amount > food_capacity:
            self.health = MAX_HEALTH
            self.family.storeFood(amount-food_capacity)
        elif amount>0:
            self.health += amount/(self.caloric_need)
            
    def individualForage(self):
        expenditure = self.current_ability*(1 - self.attributes[coll])*self.caloric_need
        self.health -= (1 - self.attributes[coll])
        food_acquired = self.location.forage(expenditure)
        self.shareWithNuclearFamily(food_acquired)
        
    def checkHunger(self):
        if self.health < 1:
            self.family.eatFoodStored(self)
        
    def marry(self, other):
        self.spouses.append(other)
        other.spouses.append(self)
        self.family += other.family
        other.family = self.family
        for child in self.children:
            other.children.append(child)
        for child in other.children:
            self.children.append(child)
            child.family = self.family
        
    def isBachelorette(self):
        return self.gender==FEMALE and not self.spouses and self.age < AGE_SENIOR and self.age >= AGE_ADULTHOOD
    
    def isNuclear(self, other):
        return set(self.parents) & set(other.parents) or self in other.parents or self in other.children
    
    def findMate(self):
        """
        STUB - find the youngest single woman
        """
        eligible_females = [f for f in self.group2 if f.isBachelorette() and not self.isNuclear(f)]
        if eligible_females:
            wife = min(eligible_females, key = lambda f: f.age)
            self.marry(wife)
        
    def doMating(self):
        if self.gender==MALE:
            if self.age >= AGE_ADULTHOOD and not self.spouses:
                self.findMate()
        else:
            if self.pregnancy>0:
                self.pregnancy += 1
                if self.pregnancy==9:
                    self.pregnancy = 0
                    ThePeople.addNewPerson(self.child, self.group2, self.location)
                    self.children.append(self.child)
                    if self.spouses: #assumes no infidelity
                        self.spouses[0].children.append(self.child)
                    self.child.family = self.family
                    #print("a birth")
            elif self.spouses and self.age <= AGE_SENIOR:
                if random()<FERT_PROB:
                    self.pregnancy = 1
                    self.child = Person.generateFromParents(self, self.spouses[0])
                    #print("a pregnancy")
    
    def step(self):
        """
        do all other behaviors besides food-getting
        """ 
        self.doMating()                        
    
    def grow(self):
        self.age += 1/12
        if self.age <= AGE_ADULTHOOD:
            self.caloric_need = 2000*30*self.age/AGE_ADULTHOOD
            if self.age >= 5:
                self.current_ability = self.ability*((self.age-5)/(AGE_ADULTHOOD-5))
            else:
                self.current_ability = 0
                
    def changeGroup(self, new_group):
        self.group2.remove(self)
        self.group2 = new_group
        new_group.append(self)

class Family(object):
    
    def __init__(self):
        self.food_stored = 0
        
    def storeFood(self, amount):
        self.food_stored = amount
        
    def eatFoodStored(self, member):
        amount = self.food_stored
        self.food_stored = 0
        member.eat(amount)
        
    def __add__(self, other):
        new_family = Family()
        new_family.food_stored = self.food_stored + other.food_stored
        return new_family
        
class Population(object):
    
    def __init__(self, population, start_location):
        self.groups = []
        group = []
        males = []
        females = []
        self.count = 0
        for i in range(population):
            person = Person.generateRandomAdult()
            self.addNewPerson(person, group, start_location)
            if person.gender==MALE:
                males.append(person)
            else:
                females.append(person)
        #pair up the males and females
        males = sorted(males, key=lambda man: man.age)
        females = sorted(females, key=lambda woman: woman.age)
        max_idx = min(len(males), len(females)) - 1
        for i in range(max_idx):
            males[i].marry(females[i])
            
        self.groups.append(group)
        self.time = 0
        
    def addNewPerson(self, person, group, location):
        group.append(person)
        TheWorld.placePersonAt(person, location)
        person.group2 = group
        self.count += 1
        if not person.patch_memory:
            person.patch_memory = [PatchMemory(location) for month in range(12)]
    
    def step(self, time):
        self.time = time
        shuffle(self.groups)
        # age all people
        for group in self.groups:
            for person in group:
                person.grow()
        
        #forage
        for group in self.groups:
            shuffle(group)
            self.group_forage(group)
            for person in group:
                person.individualForage()
                
        #rely on stored foods
        for group in self.groups:
            for person in group:
                person.checkHunger()
        
        self.cleanup()
        for group in self.groups:
            if group:
                temp = anneal(group)
                self.evaluateMove(group)
        
        # individual steps
        for group in self.groups:
            for person in group:
                person.step()
            self.group_fission(group)
        
    def cleanup(self):
        empty_groups = []
        for group in self.groups:
            dead = []
            for person in group:
                if person.health <= 0 or person.age > MAX_AGE:
                    dead.append(person)
                    
            for person in dead:
                group.remove(person)
                TheWorld.removePersonFrom(person, person.location)
                person.dies()
                self.count -= 1
            
            if not group:
                empty_groups.append(group)
                
        for group in empty_groups:
            if not group:
                self.groups.remove(group) 
                print("group died")         
            
                
    def group_forage(self, group):
        collective_fraction = 0
        for person in group:
                #first get collective action expenditures
            collective_fraction += person.attributes[coll]*person.current_ability*person.caloric_need/(2000*30)
            person.health -= person.attributes[coll]
        synergy = synergize(collective_fraction)
        #print(synergy)
        coll_expenditure = sum([synergy*person.caloric_need*person.current_ability for person in group])
        coll_food = TheWorld.patches[(0,0)].forage(coll_expenditure)
        #redistribute collectively-gained food
        portion = coll_food/len(group)
        for person in group:
            person.eat(portion)
            
    def evaluateMove(self, group, consider_current=True):
        """
        stub function - move group randomly
        """
        delta = (0,0)
        current_loc = group[0].location
        neighbors = TheWorld.neighborsOf(current_loc)
        if not consider_current:
            neighbors.remove(current_loc)
        new_loc = max(neighbors, key=lambda patch: patch.current_biomass)
        self.moveGroupTo(group, new_loc)
        
    def moveGroupTo(self, group, new_loc): 
        old_loc = group[0].location           
        for person in group:
            TheWorld.movePersonTo(person, new_loc)
            if new_loc != old_loc:
                person.family.food_stored = 0
            if new_loc.current_biomass > person.patch_memory[TheWorld.month].biomass:
                new_memory = PatchMemory(new_loc)
                person.patch_memory[TheWorld.month] = new_memory
            
    def group_fission(self, group):
        adults = [person for person in group if person.age >= AGE_ADULTHOOD]
        size = len(adults)
        try:
            max_stress = max([size/adult.attributes[grp2size] for adult in adults])
        except:
            max_stress = 0
        if max_stress > 1:  # max tolerance exceeded
            adults = sorted(adults, key = lambda adult : adult.attributes[grp2size])
            splitters = set([])  #set  
            idx = 0
            while len(splitters) <= len(adults)/2: #split in half
                idx -= 1
                splitters = splitters | set([adults[idx]]) | set(adults[idx].spouses)
            new_group = []
            for splitter in splitters:
                splitter.changeGroup(new_group)
                for child in splitter.children:
                    child.changeGroup(new_group)                
            self.groups.append(new_group)
            self.evaluateMove(new_group, consider_current=False)  
            print("group fission, new group size = %i, old_group size = %i"%(len(new_group), len(group)))            

class PatchMemory(object):
    
    def __init__(self, place):
        self.patch = place
        self.biomass = place.current_biomass
'''
Main Procedure
'''
if __name__ == '__main__':
    eff_temp = 24
    spatial_variance = 0.25
    temporal_variance = 1.0
    TheWorld = World(WORLD_DIM, World.generate, eff_temp, spatial_variance, temporal_variance)
    print("Effective temperature: %i"%eff_temp)
    print("Spatial Variance: %f"%spatial_variance)
    print("Temporal Variance: %f"%temporal_variance)
    start_location = TheWorld.patches[(0,0)]
    ThePeople = Population(INITIAL_POPULATION, start_location)   
    
    for t in range(6000):
        TheWorld.step(t)
        ThePeople.step(t)
        if t%12==0:
            locs = set([g[0].location for g in ThePeople.groups])
            dist = max([loc-start_location for loc in locs])
            print("year %i \t population %i \t groups %i \t locations %i \t distance %f"%(int(t/12), ThePeople.count, len(ThePeople.groups), len(locs), dist))
        if ThePeople.count==0:
            break
        