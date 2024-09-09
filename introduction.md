# 修改说明

## q2 random
参见代码：gui/mrc_disturb_2.py

核心思想是读取一个MRC文件，对其数据进行微扰（在原始数据大于0.1的位置加入符合正态分布的随机数），
并将结果保存在一个新的MRC文件中。

## q3 Gauss & Random
入口代码参见：gui/gd.py

```python
PROTEINS_LIST = ["/root/autodl-tmp/polnet/data/in_10A/1bxn_10A.pns"]
PROTEINS_LIST_METHOD = [1, 0]
```

- PROTEINS_LIST： 蛋白质文件对应的pns
- PROTEINS_LIST_METHOD： 采用的方法决定插入，1代表gauss，0代表random

对应调用具体行：gui/core/all_features2.py
```text
for p_id, p_file in enumerate(PROTEINS_LIST):
    ...
    net_sawlc.build_network(PROTEINS_LIST_METHOD[p_id])
```

核心逻辑实现如下：
```python
	def build_network(self, method_type=0):
		"""
		Add polymers following SAWLC model until an occupancy limit is passed
		:return:
		"""

		c_try = 0
		self._Network__pmer_fails = 0
		if self.__rots is not None:
			rot_id = self.__rot_id

		# Define the mean and standard deviation for Gaussian distribution
		mean = np.array([self._Network__voi.shape[0] / 2,
						 self._Network__voi.shape[1] / 2,
						 self._Network__voi.shape[2] / 2]) * self._Network__v_size
		std_dev = np.array([self._Network__voi.shape[0],
							self._Network__voi.shape[1],
							self._Network__voi.shape[2]]) * self._Network__v_size / 6

		# Network loop
		while (c_try <= self.__tries_pmer) and (self._Network__pl_occ < self.__occ):

			# Polymer initialization
			c_try += 1
			if self.__poly:
				p0 = np.asarray(self.__poly.GetPoint(random.randint(0, self.__poly.GetNumberOfPoints())))
			elif method_type == 0:
				p0 = np.asarray((self._Network__voi.shape[0] * self._Network__v_size * random.random(),
								 self._Network__voi.shape[1] * self._Network__v_size * random.random(),
								 self._Network__voi.shape[2] * self._Network__v_size * random.random()))
			elif method_type == 1:
				# Generate p0 using Gaussian distribution
				p0 = np.random.normal(mean, std_dev)
				p0 = np.clip(p0, [0, 0, 0],
							 np.array(self._Network__voi.shape) * self._Network__v_size)

```

## q5 细胞器插入
细胞器对应的代码入口文件参见：gui/gd.py
```python
ORGANELLE_LIST = ["/root/autodl-tmp/polnet/data/xbq/2c4c.pns"]
```

这里对应的是细胞器的pns文件，从mrc转pns，使用gui/create_macromolecule_models.ipynb 进行转换。

细胞器在插入的时候遵循的原则是首先插入，且只有一个，对应的入口函数见：gui/core/all_features2.py

```text
        # Organelle loop
        count_organelle = 0
        model_surfs, models, model_masks, model_codes = list(), list(), list(), list()
        ORGANELLE_VOI_VSIZE = VOI_VSIZE * 50
        for p_id, p_file in enumerate(ORGANELLE_LIST):
            print('\tPROCESSING FILE:', p_file)

            # Loading the protein
            organelle = MmerFile(p_file)
            # Generating the occupancy
            hold_occ = organelle.get_pmer_occ()
            if hasattr(hold_occ, '__len__'):
                hold_occ = OccGen(hold_occ).gen_occupancy()

            try:
                model = lio.load_mrc(organelle.get_mmer_svol())
            except FileNotFoundError:
                model = lio.load_mrc(ROOT_PATH + '/' + organelle.get_mmer_svol())

            model = lin_map(model, lb=0, ub=1)
            model = vol_cube(model)
            model_mask = model < organelle.get_iso()
            model[model_mask] = 0
            model_surf = pp.iso_surface(model, organelle.get_iso(), closed=False, normals=None)
            if SURF_DEC is not None:
                model_surf = pp.poly_decimate(model_surf, SURF_DEC)

            center = .5 * np.asarray(model.shape, dtype=float)
            # Monomer centering
            model_surf = pp.poly_translate(model_surf, -center)
            # Voxel resolution scaling
            model_surf = pp.poly_scale(model_surf, ORGANELLE_VOI_VSIZE)
            model_surfs.append(model_surf)
            surf_diam = pp.poly_diam(model_surf) * organelle.get_pmer_l()
            models.append(model)
            model_masks.append(model_mask)
            model_codes.append(organelle.get_mmer_id())

            # Network generation
            pol_l_generator = PGenHelixFiber()
            net_sawlc = NetOrgan(voi, ORGANELLE_VOI_VSIZE, organelle.get_pmer_l() * surf_diam, model_surf,
                                 organelle.get_pmer_l_max(),
                                 pol_l_generator, hold_occ, organelle.get_pmer_over_tol(), poly=None,
                                 svol=model < organelle.get_iso(), tries_mmer=MMER_TRIES, tries_pmer=PMER_TRIES)
            net_sawlc.build_network()

            # Density tomogram updating
            net_sawlc.insert_density_svol(model_mask, voi, ORGANELLE_VOI_VSIZE, merge='min')
            net_sawlc.insert_density_svol(model, tomo_den, ORGANELLE_VOI_VSIZE, merge='max')
            hold_lbls = np.zeros(shape=tomo_lbls.shape, dtype=np.float32)
            net_sawlc.insert_density_svol(np.invert(model_mask), hold_lbls, ORGANELLE_VOI_VSIZE, merge='max')
            tomo_lbls[hold_lbls > 0] = entity_id
            count_organelle += net_sawlc.get_num_mmers()
            print(">>>>>>>>count_organelle", count_organelle)
            og_voxels += (tomo_lbls == entity_id).sum()
            hold_vtp = net_sawlc.get_vtp()
            hold_skel_vtp = net_sawlc.get_skel()
            pp.add_label_to_poly(hold_vtp, entity_id, 'Entity', mode='both')
            pp.add_label_to_poly(hold_skel_vtp, entity_id, 'Entity', mode='both')
            pp.add_label_to_poly(hold_vtp, LBL_CP, 'Type', mode='both')
            pp.add_label_to_poly(hold_skel_vtp, LBL_CP, 'Type', mode='both')
            if poly_vtp is None:
                poly_vtp = hold_vtp
                skel_vtp = hold_skel_vtp
            else:
                poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
                skel_vtp = pp.merge_polys(skel_vtp, hold_skel_vtp)
            synth_tomo.add_network(net_sawlc, 'SAWLC', entity_id, code=organelle.get_mmer_id())
            entity_id += 1
```

- ORGANELLE_VOI_VSIZE: 对应细胞器缩放后的尺寸规模
- NetOrgan： 对应细胞器类
- build_network： 重写子类方法

具体build_network 的实现见：polnet/network.py

核心思想点：要确保对象只插入一次。
