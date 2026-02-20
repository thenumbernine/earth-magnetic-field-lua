#!/usr/bin/env luajit
--[[
The calc_b.shader takes my 10yo nvidia computer about 1 minute to link when using nMax == #wmm == 12
It will probably load a high resolution FITS file of the B field much quicker.
So let's just save one here.

I did it lazy.
- you could speed it up by using cached neighbors
- you could speed it up by multithreading it
- you could speed it up by writing it as a compute-gpu kernel (tho the nvidia linker is taking 1minute to link a glsl program made for nMax=#wmm=12. ....)
--]]
local assert = require 'ext.assert'
local cmdline = require 'ext.cmdline'(...)
local WMM = require 'earth-magnetic-field.wmm'

local W = WMM{
	cof = cmdline.cof,
	nMax = cmdline.nMax,
}
local wgs84 = W.wgs84

local Image = require 'image'

local londim = 4096
local latdim = 2048

local dphi = math.pi / londim
local dlambda = 2 * math.pi / latdim
local dheight = 1000

local BImg = Image(londim, latdim, 4, 'float')	-- Bx, By, Bz, 0
BImg:save'B.fits'
local B2Img = Image(londim, latdim, 4, 'float')	-- div B, div2D B, curl B, |curl B|
local B3Img = Image(londim, latdim, 4, 'float')	-- |B| |B|2d, 0, 0

local lastTime = os.time()
local Bptr = BImg.buffer+0
for j=0,latdim-1 do
	local v = (j + .5) / latdim
	local lat = (v * 2 - 1) * 90
	local phi = math.rad(lat)
	local cos_phi = math.cos(phi)
	for i=0,londim-1 do
		local u = (i + .5) / londim
		local lon = (u * 2 - 1) * 180
		local lambda = math.rad(lon)

		local thisTime = os.time()
		if lastTime ~= thisTime then
			lastTime = thisTime
			print(
				(100 * (i + londim * j) / (latdim * londim))
				..'%% complete'
			)
		end

		-- calcB

		local Bx, By, Bz = W:calcB(phi, lambda, 0)

		Bptr[0] = Bx	Bptr=Bptr+1
		Bptr[0] = By	Bptr=Bptr+1
		Bptr[0] = Bz	Bptr=Bptr+1
		Bptr[0] = 0		Bptr=Bptr+1
	end
end
assert.eq(Bptr, BImg.buffer + BImg.channels * BImg.width * BImg.height)


local lastTime = os.time()
local Bptr = BImg.buffer+0
local B2ptr = B2Img.buffer+0
local B3ptr = B3Img.buffer+0
for j=0,latdim-1 do
	local v = (j + .5) / latdim
	local lat = (v * 2 - 1) * 90
	local phi = math.rad(lat)
	local cos_phi = math.cos(phi)
	for i=0,londim-1 do
		local u = (i + .5) / londim
		local lon = (u * 2 - 1) * 180
		local lambda = math.rad(lon)

		local thisTime = os.time()
		if lastTime ~= thisTime then
			lastTime = thisTime
			print(
				(100 * (i + londim * j) / (latdim * londim))
				..'%% complete'
			)
		end

		local Bx, By, Bz = Bptr[0], Bptr[1], Bptr[2]

		-- calcB2

		-- [[ calculated again
		local Bx_phiR, By_phiR, Bz_phiR = W:calcB(phi + dphi, lambda, 0)
		local Bx_phiL, By_phiL, Bz_phiL = W:calcB(phi - dphi, lambda, 0)

		local Bx_lambdaR, By_lambdaR, Bz_lambdaR = W:calcB(phi, lambda + dlambda, 0)
		local Bx_lambdaL, By_lambdaL, Bz_lambdaL = W:calcB(phi, lambda - dlambda, 0)
		--]]
		--[[ use cached BImg copy
		local jR = j + 1
		local Bx_phiR, By_phiR, Bz_phiR
		if jR == latdim then
			-- across the pole?
			-- flip everything
			error"I'm lazy"
			Bx_phiR = B.buffer[0 + B.channels * (i + londim * jR)]
			By_phiR = B.buffer[1 + B.channels * (i + londim * jR)]
			Bz_phiR = B.buffer[2 + B.channels * (i + londim * jR)]
		else
			Bx_phiR = B.buffer[0 + B.channels * (i + londim * jR)]
			By_phiR = B.buffer[1 + B.channels * (i + londim * jR)]
			Bz_phiR = B.buffer[2 + B.channels * (i + londim * jR)]
		end
		--]]

		-- there is no separate altitude cache
		local Bx_heightR, By_heightR, Bz_heightR = W:calcB(phi, lambda, 0 + dheight)
		local Bx_heightL, By_heightL, Bz_heightL = W:calcB(phi, lambda, 0 - dheight)

		local dBx_dphi = (Bx_phiR - Bx_phiL) / (2 * dphi) / (wgs84.a * 1e+3 * cos_phi)
		local dBy_dphi = (By_phiR - By_phiL) / (2 * dphi) / (wgs84.a * 1e+3 * cos_phi)
		local dBz_dphi = (Bz_phiR - Bz_phiL) / (2 * dphi) / (wgs84.a * 1e+3 * cos_phi)

		local dBx_dlambda = (Bx_lambdaR - Bx_lambdaL) / (2 * dlambda) / (wgs84.a * 1e+3)
		local dBy_dlambda = (By_lambdaR - By_lambdaL) / (2 * dlambda) / (wgs84.a * 1e+3)
		local dBz_dlambda = (Bz_lambdaR - Bz_lambdaL) / (2 * dlambda) / (wgs84.a * 1e+3)

		local dBx_dheight = (Bx_heightR - Bx_heightL) / (2 * dheight) / (wgs84.a * 1e+3)
		local dBy_dheight = (By_heightR - By_heightL) / (2 * dheight) / (wgs84.a * 1e+3)
		local dBz_dheight = (Bz_heightR - Bz_heightL) / (2 * dheight) / (1e+3)

		local div2D_B = dBx_dphi + dBy_dlambda
		local div_B = div2D_B + dBz_dheight

		local curl_B_x = dBz_dlambda - dBy_dheight
		local curl_B_y = dBx_dheight - dBz_dphi
		local curl_B_z = dBy_dphi - dBx_dlambda

		local len_curl_B = math.sqrt(curl_B_x^2 + curl_B_y^2 + curl_B_z^2)

		B2ptr[0] = div_B		B2ptr=B2ptr+1
		B2ptr[0] = div2D_B		B2ptr=B2ptr+1
		B2ptr[0] = curl_B_z		B2ptr=B2ptr+1
		B2ptr[0] = len_curl_B	B2ptr=B2ptr+1


		-- calcB3


		local len_B2 = math.sqrt(Bx^2 + By^2)
		local len_B3 = math.sqrt(len_B2^2 + Bz^2)

		B3ptr[0] = len_B2		B3ptr=B3ptr+1
		B3ptr[0] = len_B3		B3ptr=B3ptr+1
		B3ptr[0] = 0			B3ptr=B3ptr+1
		B3ptr[0] = 0			B3ptr=B3ptr+1


		Bptr=Bptr+4
	end
end
assert.eq(Bptr, BImg.buffer + BImg.channels * BImg.width * BImg.height)
assert.eq(B2ptr, B2Img.buffer + B2Img.channels * B2Img.width * B2Img.height)
assert.eq(B3ptr, B3Img.buffer + B3Img.channels * B3Img.width * B3Img.height)

BImg:save'B.fits'
B2Img:save'B2.fits'
B3Img:save'B3.fits'
