
from bs4 import BeautifulSoup
import requests
import urllib2

num=int(80)

url = 'http://www.astro.umd.edu/~jph/KURUCZ_CONTINUUM/'
ext = '-7464A'
page = requests.get(url).text
soup = BeautifulSoup(page, 'html.parser')
address=[url  + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]
filename=[node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

count= len(filename)
print "tottal number of files should be downlouded: ", count-num

print "count: ", count
input("Enter a number ")


for i in range(count-num):
    print "Downloading the file number: ", i+num, "total: ", count
    print "its addess:" , address[i+num]
    print "its name: ", str(filename[i+num])
    f = open(str(filename[i+num]),'wb')
    f.write(urllib2.urlopen(address[i+num]).read())
    f.close()
    print "************************************************************"


