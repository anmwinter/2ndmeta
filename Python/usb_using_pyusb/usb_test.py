##hopefully an attempt to get data off a datalogger co2meter.com model ZGm053U
##

##You will need to know your Vendor ID and Product ID 

import usb.core
import usb.util
import sys
import logging

##To debug
## env PYUSB_DEBUG_LEVEL=debug python usb_test.py
interface = 0

dev = usb.core.find()

# find the USB device
dev = usb.core.find(idVendor = 0x04d9, idProduct = 0xa052)

if dev is None:
    print 'Device not found'+'\n'
else:
    print 'USB CO2 and Temperature Sensor Found'

for cfg in dev:
    sys.stdout.write(str(cfg.bMaxPower) + " mA" +'\n')
    sys.stdout.write(str(cfg.bDescriptorType) +'\n')
    sys.stdout.write(str(cfg.bmAttributes) +'\n')

if dev.is_kernel_driver_active(interface) is True:
        print "but we need to detach kernel driver"
        dev.detach_kernel_driver(interface)
        print "claiming device"
        usb.util.claim_interface(dev, interface)
else: print 'Not claimed!!'

dev_mem = usb.core.find(idVendor = 0x05ac, idProduct = 0x8404)

if dev_mem is None:
    print 'Device not found'
else:
    print 'Internal Memory Found On Sensor Found'

for cfg in dev_mem:
    sys.stdout.write(str(cfg.bMaxPower) + " mA" +'\n')
    sys.stdout.write(str(cfg.bDescriptorType) +'\n')
    sys.stdout.write(str(cfg.bmAttributes) +'\n')



