import urllib
import pywikibot
from pywikibot import pagegenerators
import requests
from PIL import Image
import imagehash
import io


import math
import imagehash
import scipy.fftpack
import numpy
from PIL import Image, ImageFilter
import numpy as np
from PIL import Image

# Define the precision bits equivalent to the C++ #define
PRECISION_BITS = (32 - 8 - 2)

# Create the lookup table for values from -640 to 639
clip8_lookups = [0] * 640 + list(range(256)) + [255] * 384

def imaging_new_dirty(mode, xsize, ysize):
    # Create a new empty image
    return Image.new(mode, (xsize, ysize))

def sinc_filter(x):
    """Compute the sinc filter value."""
    if x == 0.0:
        return 1.0
    x = x * math.pi
    return math.sin(x) / x

def lanczos_filter(x):
    """Compute the Lanczos filter value, which is a truncated sinc function."""
    if -3.0 <= x < 3.0:
        return sinc_filter(x) * sinc_filter(x / 3)
    return 0.0

def clip8(in_value):
    """
    Clips the input integer to the range [0, 255] using a lookup table.
    
    Parameters:
        in_value (int): The input integer value to clip.
        
    Returns:
        int: The clipped value in the range [0, 255].
    """
    # Compute the index by shifting the input value
    index = (int(in_value) >> PRECISION_BITS) + 640
#    index = round(in_value) >> PRECISION_BITS + 640
    # Return the corresponding value from the lookup table
    print(f"FOO3: {in_value} {index} {clip8_lookups[index]}")
    return clip8_lookups[index]


def normalize_coeffs_8bpc(out_size, ksize, prekk):
    """
    Normalize coefficients similar to the C++ implementation.

    Parameters:
    - out_size: int, the output size.
    - ksize: int, the kernel size.
    - prekk: numpy array of double (float64), the precomputed coefficients.

    Returns:
    - kk: numpy array of int32, the normalized coefficients.
    """
    # Allocate an array of int32 for the normalized coefficients
    kk = np.zeros_like(prekk, dtype=np.int32)

    # Perform normalization
    for x in range(out_size * ksize):
#        print(f"PREKK_FOO:{x} {prekk[x]}")
        if prekk[x] < 0:
            kk[x] = int(-0.5 + prekk[x] * (1 << PRECISION_BITS))
#            kk[x] = round(-0.5 + prekk[x] * (1 << PRECISION_BITS))
        else:
            kk[x] = int(0.5 + prekk[x] * (1 << PRECISION_BITS))
#            kk[x] = round(0.5 + prekk[x] * (1 << PRECISION_BITS))

    return kk


def precompute_coeffs(inSize, in0, in1, outSize, filterp):
    support, scale, filterscale = 0.0, 0.0, 0.0
    center, ww, ss = 0.0, 0.0, 0.0
    xx, x, ksize, xmin, xmax = 0, 0, 0, 0, 0

    # prepare for horizontal stretch
    filterscale = scale = (in1 - in0) / outSize
    if filterscale < 1.0:
        filterscale = 1.0

    # determine support size (length of resampling filter)
    support = filterp['support'] * filterscale

    # maximum number of coeffs
    ksize = math.ceil(support) * 2 + 1

    # check for overflow
    if outSize > (2**31 - 1) / (ksize * 8):  # 8 is sizeof(double) in most systems
        raise MemoryError("Integer overflow")

    # coefficient buffer
    kk = np.zeros(outSize * ksize, dtype=np.float64)
    bounds = np.zeros(outSize * 2, dtype=np.int32)

    for xx in range(outSize):
        center = in0 + (xx + 0.5) * scale
        ww = 0.0
        ss = 1.0 / filterscale

        # Calculate xmin and xmax
        xmin = max(int(center - support + 0.5), 0)
        xmax = min(int(center + support + 0.5), inSize)
#        xmin = max(round(center - support + 0.5), 0)
#        xmax = min(round(center + support + 0.5), inSize)
        xmax -= xmin

        # Compute weights for the current output pixel
        offset = xx * ksize
        for x in range(xmax):
            w = filterp['filter']((x + xmin - center + 0.5) * ss)
            kk[offset + x] = w
            ww += w

        for x in range(xmax):
            if ww != 0.0:
                kk[offset + x] /= ww

        # Remaining values should stay as zero if they are used despite of xmax
        for i in range(x, ksize):
            kk[offset + i] = 0

        bounds[xx * 2] = xmin
        bounds[xx * 2 + 1] = xmax


    print(f'precompute_coeffs: {ksize}\n{bounds}\n{kk}')
    return ksize, bounds, kk

def precompute_coeffs_old(in_size, in0, in1, out_size, filterp):
    """
    Precompute the coefficients for resampling.
    
    Parameters:
    - in_size: int, size of the input dimension
    - in0: float, start of the input range
    - in1: float, end of the input range
    - out_size: int, size of the output dimension
    - filterp: dict, containing 'filter' function and 'support' value

    Returns:
    - ksize: int, number of coefficients per output pixel
    - bounds: numpy array (1D), the bounds for each output pixel
    - kk: numpy array (1D), the precomputed coefficients
    """
    # Prepare for horizontal stretch
    print(f'IN: {in1}\t{in0}\t{out_size}')

    scale = (in1 - in0) / out_size
    filterscale = max(1.0, scale)
    
    # Determine support size (length of resampling filter)
    support = filterp['support'] * filterscale
    print(filterp['support'])
    print(filterscale)
    
    # Maximum number of coefficients
    ksize = int(math.ceil(support) * 2 + 1)

    print(f'{in_size}\t{in0}\t{in1}\t{out_size}\t{support}\t{ksize}\n')
    
    # Initialize coefficient buffer and bounds
    kk = np.zeros(out_size * ksize, dtype=np.float64)
    bounds = np.zeros(out_size * 2, dtype=np.int32)

    for xx in range(out_size):
        center = in0 + (xx + 0.5) * scale
        ww = 0.0
        ss = 1.0 / filterscale
        
        # Calculate xmin and xmax
        xmin = max(int(center - support + 0.5), 0)
        xmax = min(int(center + support + 0.5), in_size)
        xmax -= xmin
        
        # Compute weights for the current output pixel
        offset = xx * ksize
        for x in range(xmax):
            w = filterp['filter']((x + xmin - center + 0.5) * ss)
            kk[offset + x] = w
            ww += w
        
        # Normalize the weights
        if ww != 0.0:
            kk[offset:offset + xmax] /= ww
        
        # Fill the remaining coefficients with zeros if not used
        kk[offset + xmax:offset + ksize] = 0
        
        # Store the bounds in 1D array
        bounds[xx * 2] = xmin
        bounds[xx * 2 + 1] = xmax
    
    return ksize, bounds, kk


def imaging_resample_horizontal_8bpc(im_out, im_in, offset, ksize, bounds, prekk):

#    def make_uint32(a, b, c, d):
#        return (a << 24) | (b << 16) | (c << 8) | d

    im_out_np = np.array(im_out)  # Convert PIL image to NumPy array
    im_in_np = np.array(im_in)  # Convert PIL image to NumPy array
    im_in_bytes = im_in.tobytes()

    # Normalize coefficients (assuming this function is defined elsewhere)
    kk = normalize_coeffs_8bpc(im_out_np.shape[1], ksize, prekk)

    ysize=im_out_np.shape[0]
    xsize=im_out_np.shape[1]

#    print(f'{ysize}\t{xsize}')
    if len(im_in_np.shape) == 2:  # Assuming this is equivalent to imIn->image8
        for yy in range(ysize):
            for xx in range(xsize):
                xmin = bounds[xx * 2]
                xmax = bounds[xx * 2 + 1]
#                print(f'XMIN_XMAX: {xmin}\t{xmax}\n')
#                k = kk[xx * ksize : (xx + 1) * ksize]
                pos_k = xx * ksize
                ss0 = 1 << (PRECISION_BITS - 1)
                for x in range(xmax):
                    ss0 += im_in_np[yy + offset, x + xmin] * kk[pos_k + x]
                #    print(f'FOO\t{yy}\t{xx}\t{x}\t{im_in_np[yy + offset, x + xmin] * kk[pos_k + x]}\t{kk[pos_k + x]}\t{im_in_np[yy + offset, x + xmin]}, {im_in_bytes[yy + offset + x + xmin]}')
#                print(f"SS0 PY: {ss0}")
#                print("----")
                im_out_np[yy, xx] = clip8(ss0)
#                print(im_out_np[yy, xx]);

    return Image.fromarray(im_out_np.astype('uint8'), 'L')

def imaging_resample_vertical_8bpc(im_out, im_in, offset, ksize, bounds, prekk):
    """Performs vertical resampling on an 8-bit per channel image."""
    im_out_np = np.array(im_out)  # Convert PIL image to NumPy array
    im_in_np = np.array(im_in)  # Convert PIL image to NumPy array
    ysize, xsize = im_out_np.shape
    kk = normalize_coeffs_8bpc(ysize, ksize, prekk)

    print(f'{ysize}\t{xsize}')
    for yy in range(ysize):
        print("===");
        ymin = bounds[yy * 2 + 0]
        ymax = bounds[yy * 2 + 1]
#        k = kk[yy * ksize:(yy + 1) * ksize]
        print(f'FOO YMIN_YMAX: {ymin} {ymax}')
        for xx in range(xsize):
            pos_k = yy*ksize
            ss0 = 1 << (PRECISION_BITS - 1)

            for y in range(ymax):
                ss0 += im_in_np[y + ymin, xx] * kk[pos_k + y]
            print(f'FOO SS0: {ss0}')
            im_out_np[yy, xx] = clip8(ss0)
            print(f'FOO2: {im_out_np[yy, xx]}')


    return Image.fromarray(im_out_np.astype('uint8'), 'L')


def imaging_resample_inner(im_in, xsize, ysize, filterp, box, resample_horizontal, resample_vertical):
    """
    Perform a two-pass resize operation with horizontal and vertical passes.
    """
    need_horizontal = xsize != im_in.width or box[0] != 0 or box[2] != xsize
    need_vertical = ysize != im_in.height or box[1] != 0 or box[3] != ysize


    print(f'{xsize}\t{ysize}\n')

    # Horizontal coefficients
    ksize_horiz, bounds_horiz, kk_horiz = precompute_coeffs(
        im_in.width, box[0], box[2], xsize, filterp
    )
    if not ksize_horiz:
        return None

    # Vertical coefficients
    ksize_vert, bounds_vert, kk_vert = precompute_coeffs(
        im_in.height, box[1], box[3], ysize, filterp
    )
    if not ksize_vert:
        return None

    ybox_first = bounds_vert[0]
    ybox_last = bounds_vert[ysize * 2 - 2] + bounds_vert[ysize * 2 - 1]

    print("-------\n")
    print(f'{ybox_first}\t{ybox_last}')

    im_temp = None
    im_out = None

    # Horizontal pass
    if need_horizontal:
        for i in range(ysize):
            bounds_vert[i * 2] -= ybox_first

        im_temp = imaging_new_dirty(im_in.mode, xsize, ybox_last - ybox_first)
        if im_temp:
            im_temp=resample_horizontal(im_temp, im_in, ybox_first, ksize_horiz, bounds_horiz, kk_horiz)
        if not im_temp:
            return None
        im_out = im_temp
        print("Horiz pass")

    # Vertical pass
    if need_vertical:
        im_out = imaging_new_dirty(im_in.mode, xsize, ysize)
        if im_out:
            im_out = resample_vertical(im_out, im_temp or im_in, 0, ksize_vert, bounds_vert, kk_vert)
        if not im_out:
            return None
        print("Vertical pass")


    # If neither horizontal nor vertical passes were needed, copy the image
    if not im_out:
        im_out = im_in.copy()

    return im_out


def lanczos_resample(image, output_size, a=3):
    """Resize image using Lanczos filter, mimicking Pillow's method."""

    box=[]
    box.append(0)      
    box.append(0)      
    box.append(image.size[0])      
    box.append(image.size[1])
    (xsize, ysize) = output_size
    filterp =  {
        'filter': lanczos_filter,
        'support': 3.0
        }

    return imaging_resample_inner(
        image, xsize, ysize, filterp, box, imaging_resample_horizontal_8bpc, imaging_resample_vertical_8bpc);
 



#    image_np = np.array(image, dtype=np.int8)
    
#    ksize, bounds, coeffs = precompute_coeffs(
#        in_size=image_np.shape[1],
#        in0=0.0,
#        in1=image_np.shape[1],
#        out_size=output_size[1],
#        a=a
#    )
    
    resized_image = apply_resampling(image_np, output_size, ksize, bounds, coeffs)
    return Image.fromarray(np.clip(resized_image, 0, 255).astype(np.uint8))

# Example usage
#image = Image.open("path_to_your_image.jpg")
#output_size = (new_height, new_width)
#resized_image = lanczos_resample(image, output_size)
#resized_image.save("resized_image.jpg")
# -------------

def custom_grayscale(image):
    # Load image data
    pixels = list(image.getdata())
    width, height = image.size

    # Grayscale conversion using the standard luminance formula
    grayscale_pixels = []
    for pixel in pixels:
        # If the image is RGB, we expect a tuple of (R, G, B)
        if len(pixel) == 3:
            r, g, b = pixel
        # If the image has an alpha channel (RGBA), we ignore the alpha
        elif len(pixel) == 4:
            r, g, b, _ = pixel
        
        # Compute the grayscale value using the luminance formula
        grayscale_value = round(0.299 * r + 0.587 * g + 0.114 * b)
        
        # Store the grayscale value as a tuple (R=G=B=grayscale_value)
        grayscale_pixels.append(grayscale_value)
    
    # Create a new grayscale image
    grayscale_image = Image.new('L', (width, height))
    grayscale_image.putdata(grayscale_pixels)

    return grayscale_image


def phash1(image, hash_size=8, highfreq_factor=4):
    img_size = hash_size * highfreq_factor
    image_grayscale = image.convert('L')
    image_resized= image_grayscale.resize((img_size, img_size), Image.Resampling.LANCZOS)

    # Convert the grayscale image to a byte array
    grayscale_bytes = image_grayscale.tobytes()

    # Print the first 64 bytes
    hex_values = grayscale_bytes[:64].hex()
#    print(hex_values)

    resized_bytes = image_resized.tobytes()
    hex_values = resized_bytes[:256].hex()
    print("")
    print(hex_values)

    pixels = numpy.asarray(image_resized)
    dct = scipy.fftpack.dct(scipy.fftpack.dct(pixels, axis=0), axis=1)
    dctlowfreq = dct[:hash_size, :hash_size]
    med = numpy.median(dctlowfreq)
    diff = dctlowfreq > med
    return imagehash.ImageHash(diff)

def phash2(image, hash_size=8, highfreq_factor=4):
    img_size = hash_size * highfreq_factor
    image_grayscale = custom_grayscale(image)
#    image_grayscale = image.convert('L')

#    image_resized= image_grayscale.resize((img_size, img_size), Image.Resampling.LANCZOS)
    image_resized = lanczos_resample(image_grayscale, (img_size, img_size))

    # Convert the grayscale image to a byte array
    grayscale_bytes = image_grayscale.tobytes()

    # Print the first 64 bytes
    hex_values = grayscale_bytes[:64].hex()
#    print(hex_values)

    resized_bytes = image_resized.tobytes()
    hex_values = resized_bytes[:256].hex()
    print("")
    print(hex_values)


    pixels = numpy.asarray(image_resized)
    dct = scipy.fftpack.dct(scipy.fftpack.dct(pixels, axis=0), axis=1)
    dctlowfreq = dct[:hash_size, :hash_size]
    med = numpy.median(dctlowfreq)
    diff = dctlowfreq > med
    return imagehash.ImageHash(diff)

image = Image.open('test.jpg')
t1=phash1(image)
t2=phash2(image)


print(t1)
print(t2)



site = pywikibot.Site('commons', 'commons')
generator = pagegenerators.RandomPageGenerator(total=50, site=site, namespaces=[6])

images = []
for page in generator:
    if page.title().lower().endswith('.jpg'):
        images.append(page)
    if len(images) >= 10:
        break

# Process each image: get URL, download the image, and calculate pHash
for image in images:
    image_info = image.latest_file_info
    image_url = image_info.url  # Get the URL of the image file
    print(f"Image title: {image.title()}")
    print(f"Image URL: {image_url}")

    # Download the image
    image_file = Image.open(urllib.request.urlopen(image_url))

    # Calculate perceptual hash (pHash)
    p1 = imagehash.phash(image_file)
    p2 = phash2(image_file)
    if p1 != p2: 
        print(f"Image title: {image.title()}")
        print(f"Image URL: {image_url}")
        print(f"pHash1: {p1}\n")
        print(f"pHash2: {p2}\n")
        print("VIRHE")
        exit(1)
  
print("OK")
