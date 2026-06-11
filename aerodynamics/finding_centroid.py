import cv2
import numpy as np

# img = cv2.imread('aerodynamics\Side_view.jpeg', cv2.IMREAD_GRAYSCALE)
# _, thresh = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY_INV)
# M = cv2.moments(thresh)
# cx = int(M['m10'] / M['m00'])
# cy = int(M['m01'] / M['m00'])

import cv2
import numpy as np

img = cv2.imread('aerodynamics/Side_view.jpeg')
print(img.shape)  # (height, width, channels)
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

_, thresh = cv2.threshold(gray, 200, 255, cv2.THRESH_BINARY_INV)

contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
largest = max(contours, key=cv2.contourArea)

mask = np.zeros_like(gray)
cv2.drawContours(mask, [largest], -1, 255, thickness=cv2.FILLED)

M = cv2.moments(mask)
cx = int(M['m10'] / M['m00'])
cy = int(M['m01'] / M['m00'])

print(f"Centroid: x={cx}, y={cy}")

# Draw and save result
result = img.copy()
cv2.circle(result, (cx, cy), 12, (0, 0, 255), -1)
cv2.line(result, (cx-25, cy), (cx+25, cy), (255,255,255), 2)
cv2.line(result, (cx, cy-25), (cx, cy+25), (255,255,255), 2)
cv2.imwrite('centroid_result.png', result)

height_minus_lg = img.shape[0] / img.shape[1] * 11
centroidx = cx / img.shape[1] * 11
centroidy = cy / img.shape[0] * height_minus_lg
print(f'centroid: (x,z) = {[centroidx, centroidy]}')
dist_cg_to_centroid = np.sqrt((centroidx - 17.921 * 0.3048)**2 + (7.906 * 0.3048 - centroidy)**2)
print(f'Distance cg to centroid = {dist_cg_to_centroid}')


# Class II  computed c.g.:  x=17.921 ft,  y=0.000 ft,  z=5.906 ft